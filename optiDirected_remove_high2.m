%A software package for optimizing synchronization of coupled oscillators with high-order networks
%(c) 2021 Ying Tang, Dinghua Shi, Linyuan Lv
%All rights reserved. 
%This MATLAB code package optimizes network topology for synchronization of coupled oscillators 
%with high-order interactions. The current focus is the system with Kuramoto-type coupling function 
%for identical oscillators, the second-order interactions (triangle). The optimization is realized by 
%minimizing the eigenratios or the spread of eigenvalues for the generalized Laplacian matrices. For the undirected network, 
%we rewire the triangle interactions and use simulated annealing to optimize the network synchronizability. 
%For the directed network, we selectively remove directional triangle interactions to optimize synchronizability, 
%and investigate asymmetry for the optimized directed network.

%A detailed description on the scripts is in README file. 
%Contact: Ying Tang, jamestang23@gmail.com

%% Initilize parameters
clear;
addpath([pwd,'\octave-networks-toolbox-master']);
figurefolder=pwd;
figureSubfolder=[figurefolder,'\DirectedNetworkAsymmetry2'];
mkdir(figureSubfolder);
InitialNet=1;%the two types of initialization as described in text: chose 1 or 2.

Run =1e3;% the number of independent simulated annealings to be performed
NumIncreaseTarget=1e2;% the number of simulated annealings before the target eigenratio increased by 10%, because the target eigenratio can not be too small as it may not be reached due to the sparsity for large networks

%% Start multiple numerical runs
for realization=1:30%1:10%10:-1:6%1:10 % the number of numerical replicates
tic;
filename=[figureSubfolder,'\RemoveScan_',num2str(realization),'_InitialNet_',num2str(InitialNet),'.mat'];%Scan has target 1.8/1.4, Scan2 has target value 1


for NumberNode=[6, 10, 20, 50, 80, 100]% Scan different network sizes
     TriangleNum=NumberNode*2;
    disp(['Network size ', num2str(NumberNode)]);toc;
    ScanDelete=min(round(NumberNode*0.9),NumberNode-1)*InitialNet; % at each iteration between 1 to ScanDelete triangles will be randomly removed
    clear AdjacencyTensor;

%% Input network structure
if InitialNet==1% First type of initialized network
for j=1:NumberNode%NumberNode
    id=mod([1,2,3]+(j-1),NumberNode)+1;
    P = perms(id); %symmetric tensor for each triangle
    for k=1:size(P,1)
    AdjacencyTensor(P(k,1),P(k,2),P(k,3))=1;
    end
end
%Add randomly to have TriangleNum new triangles
[AdjConverted,triangle_initial]=ConvertTriangleToAdjacency(AdjacencyTensor);
TriangleNumCurrent=size(triangle_initial,1);
while TriangleNumCurrent<TriangleNum
    id=randperm(NumberNode,3);
     P = perms(id); %symmetric tensor for each triangle
    for k=1:size(P,1)
    AdjacencyTensor(P(k,1),P(k,2),P(k,3))=1;
    end  
    [AdjConverted,triangle_initial]=ConvertTriangleToAdjacency(AdjacencyTensor);
    TriangleNumCurrent=size(triangle_initial,1);
end
end
 

%% Calculate Laplacian and eigenvalues for high-order
% M is the number of directional links in the original network
idx = find(AdjacencyTensor == 1);
M = size(idx,1);%ddd
% construct Laplacian matrix
lap=Laplacian2(AdjacencyTensor);
lapStart=lap;
% find eigen-ratio of the original network
eigv = sort(eig(lap),'ComparisonMethod','real');
% find the normalized spread of the eigenvalues
lambdabar=sum(eigv(2:end))/(NumberNode-1);%/(NumberNode-2);
d=length(find(AdjacencyTensor == 1))/NumberNode/(NumberNode-1);
SpreadEigenvalues=1/d^2/(NumberNode-1)/(NumberNode-2)*sum(abs(eigv(2:end)-lambdabar).^2);

r_original = SpreadEigenvalues;%eigv(end)/eigv(2);
% target value of the eigen-ratio:  too small makes convergence harder
r_targetFinal = 0.1;%1;%
r_target = r_targetFinal;%


%% Start optimization

% remove records the number of links removed before the target eigen-ratio is reached at each run
remove = zeros(1,Run);
% inverse temperature used in simulated annealings (different for each run)
Beta = linspace(50,500,Run);
r_optimal=Inf;
OptimalNetwork=AdjacencyTensor;

for i = 1:Run
    %r_target = r_targetFinal;%1;%
    if mod(i,1e3)==0
        disp(i);
    end
    adj = AdjacencyTensor;
    r_new = r_original;
    r_current = r_original;
     Metropolis=0;
     beta = Beta(i);
    loopcnt = 0;
    % start of simulated annealing
    while r_new > r_target || Metropolis==0
        aij = adj;
        
        idx = find(adj == 1);
        num_rem = randi(min(ScanDelete,size(idx,1)));     
        del_link = idx(randperm(size(idx,1),num_rem));     
        aij(del_link) = 0;
        
        % find the eigen-ratio of the new network
        %module Laplacian matrix:
       lap=Laplacian2(aij);
        eigv = sort(eig(lap),'ComparisonMethod','real');
        % if the smallest nontrival eigenvalue becomes zero, exit loop
        if abs(eigv(2))<1e-5
           remove(i) = Inf;
           break
        end
        %%r_new = real(eigv(end))/real(eigv(2));
        % find the normalized spread of the eigenvalues
        %disp(eigv);
        %lambdabar=sum(eigv(2:end))/(NumberNode-1);
        %d=length(find(aij == 1))/NumberNode/(NumberNode-1);
        %SpreadEigenvalues=1/d^2/(NumberNode-1)*sum(abs(eigv(2:end)-lambdabar).^2);
        lambdabar=sum(eigv(2:end))/(NumberNode-1);%/(NumberNode-2);
        d=length(find(aij == 1))/NumberNode/(NumberNode-1);
        SpreadEigenvalues=1/d^2/(NumberNode-1)/(NumberNode-2)*sum(abs(eigv(2:end)-lambdabar).^2);

        r_new =SpreadEigenvalues;
        
        % energy is negative if the eigen-ratio is improved by removal
        energy = r_new/min(r_original,r_current)-1;
        %disp(r_target);
        %disp([SpreadEigenvalues, r_original,r_current,energy]);
       % ddd
        % Metropolis criterion
        %Metropolis=0;
        if rand < exp(-beta*energy)
            adj = aij;
            r_current = r_new;
            remove(i) = remove(i)+num_rem;
            Metropolis=1;
            % if more than M-NumberNode removals have been performed, abort this run
            if remove(i) > M-NumberNode
                remove(i) = Inf;
                break
            end
        end
        
        
        loopcnt = loopcnt + 1;%disp(loopcnt);
        if mod(loopcnt,NumIncreaseTarget)==0% increase r_target by 10% to ensure finding optimal network
            loopcnt=0;
            r_target=abs(r_target*1.1);
        end
        
        
    end
    
    
    if Run==1
        OptimalNetwork=adj;
         r_target=r_current;%my code
        r_optimal=r_current;
    end
    
    % if the target eigen-ratio is reached
    if r_current <= min(r_target,r_optimal)
         OptimalNetwork=adj;
         r_target=r_current;%my code
        r_optimal=r_current;
    end
end

OptimalLap=real(Laplacian2(OptimalNetwork));
Summary.EigenRatio(NumberNode,1)=r_original;
Summary.EigenRatio(NumberNode,2)=r_optimal;
Summary.AdjTensorInitial{NumberNode}=AdjacencyTensor;
Summary.AdjTensorOptimal{NumberNode}=OptimalNetwork;

[Summary.AsymmetryIndex1{NumberNode},Summary.AsymmetryIndex2{NumberNode}]=AsymmetryIndex(OptimalNetwork);

end
save(filename,'Summary');
end
%%
function [AsymmetryIndex1, AsymmetryIndex2]=AsymmetryIndex(OptimalNetwork)

AsymmetryIndex1=sum(sum(abs(OptimalNetwork-permute(OptimalNetwork,[1,3,2])),3),2);
Temp1=(permute(OptimalNetwork,[2,1,3])+permute(OptimalNetwork,[2,3,1]));
Temp2=(permute(OptimalNetwork,[3,1,2])+permute(OptimalNetwork,[3,2,1]));
Temp3=(permute(OptimalNetwork,[1,2,3])+permute(OptimalNetwork,[1,3,2]));

AsymmetryIndex2=sum(sum(abs(Temp1+Temp2-2*Temp3),3),2);
end

function lap=Laplacian1(AdjacencyTensor)
rowsum = sum(AdjacencyTensor,2);
lap = diag(rowsum) - AdjacencyTensor;
end