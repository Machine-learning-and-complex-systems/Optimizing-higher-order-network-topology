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
figureSubfolder=[figurefolder,'\ScanSecondOrder5'];
mkdir(figureSubfolder);
InitialNet=1;%the two types of initialization as described in text: chose 1 or 2.

Run =1e3;% the number of independent simulated annealings to be performed
NumIncreaseTarget=1e2;% the number of simulated annealings before the target eigenratio increased by 10%, because the target eigenratio can not be too small as it may not be reached due to the sparsity for large networks


%% Start multiple numerical runs
for realization=1%:10% the number of numerical replicates

filename=[figureSubfolder,'\Rewire2Scan4_',num2str(realization),'_InitialNet_',num2str(InitialNet),'.mat'];%Scan has target 1.8/1.4, Scan2 has target value 1

tic;
for NumberNode=[6:10]% Scan different network sizes
    
    TriangleNum=NumberNode*2;%Fix the number of triangles
    
    disp(['Network size ', num2str(NumberNode)]);toc;
    ScanDelete=min(round(NumberNode*0.2),NumberNode-1);% at each iteration between 1 to ScanDelete triangles will be randomly removed
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

%%fully connected network
% for j=1:NumberNode
%     for jj=1:NumberNode
%         for jjj=1:NumberNode
%             if j~=jj && jj~=jjj && j~=jjj
%                 vtg(j,jj,jjj)=1;
%             end
%         end
%     end
% end


%% Calculate Laplacian and eigenvalues for high-order network
% M is the number of directional links in the original network
idx0 = find(AdjacencyTensor == 1);
M = size(idx0,1);
% NumLinks is the number of links in the corresponding first-order network
[AdjConverted,triangle]=ConvertTriangleToAdjacency(AdjacencyTensor);
idx1 = find(AdjConverted == 1);
NumLinks=size(idx1,1)/2;
% construct Laplacian matrix
lap=Laplacian2(AdjacencyTensor);
% find eigen-ratio of the original network
eigv = sort(eig(lap));
r_original = eigv(end)/eigv(2);
% target value of the eigen-ratio:  too small makes convergence harder
r_targetFinal=1;%r_target = 1;%1.8;%r_original/2 + 0.5;
r_target = r_targetFinal;%1;%

%% Start optimization

% inverse temperature used in simulated annealings (different for each run)
Beta = linspace(1,500,Run);
% remove records the number of links removed before the target eigen-ratio is reached at each run
remove = zeros(1,Run);

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
    % start of simulated annealing
    loopcnt = 0;
    num_rewire = randi(min(ScanDelete,size(idx0,1)));   %randi(4);
        
    while r_new > r_target || Metropolis==0
        aij = adj;
        
        % at each iteration between 1 to ScanDelete triangles will be randomly rewired
        aaa=[];
        idx = find(aij == 1);
        [i1, i2, i3] = ind2sub(size(aij), idx);
        aaa(:,1)=i1;aaa(:,2)=i2;aaa(:,3)=i3;
        aaa=unique(sort(aaa,2),'rows');
        %disp(aaa);dddd
        del_triangle = randsample(length(aaa),num_rewire);%randperm(length(idx),num_rewire);%allowing to repeat deleting the same triangle    
        %disp(del_triangle);
        del_index=[aaa(del_triangle,1), aaa(del_triangle,2), aaa(del_triangle,3)];
        for kk=1:size(del_index,1)
            del_links=perms(del_index(kk,:));
            for ll=1:size(del_links,1)
            aij(del_links(ll,1),del_links(ll,2),del_links(ll,3)) = 0;
            end
        end
        
        %add triangle back to the adjacency tensor.
         addtriangle=1;bbb=[];
        while addtriangle==1
        ind = find(aij == 0);
        [ii1, ii2, ii3] = ind2sub(size(aij), ind);
        bbb(:,1)=ii1;bbb(:,2)=ii2;bbb(:,3)=ii3;
        bbb=unique(sort(bbb,2),'rows');
        bbb(any(diff(sort(bbb,2),[],2)==0,2),:)=[];
        add_triangle = randsample(length(bbb),num_rewire);
        %add_triangle = randsample(length(idx),num_rewire);%randperm(length(ind),num_rewire);%allowing to repeat deleting the same triangle      
        add_index=[bbb(add_triangle,1), bbb(add_triangle,2), bbb(add_triangle,3)];
        add_index = unique(sort(add_index,2),'rows');
        checkIdentical=diff(add_index,1,2);
        for kk=1:size(add_index,1)
            if ~any(checkIdentical(kk,:)==0)        
            add_links=perms(add_index(kk,:));
            for ll=1:size(add_links,1)
            aij(add_links(ll,1),add_links(ll,2),add_links(ll,3)) = 1;
            end
             addtriangle=0;
            end
        end
        end
              
        % find the eigen-ratio of the new network
        %module Laplacian matrix:
       lap=Laplacian2(aij);
        eigv = sort(eig(lap),'ComparisonMethod','real');
        % if the smallest nontrival eigenvalue becomes zero, exit loop
        if abs(eigv(2))<1e-5
           break
        end
        r_new = real(eigv(end))/real(eigv(2));
        % energy is negative if the eigen-ratio is improved by removal
        energy = r_new/min(r_original,r_current)-1;
       % disp(r_target);
        %disp([r_original,r_current,energy]);
        % Metropolis criterion
        %Metropolis=0;
        if rand < exp(-beta*energy)
            adj = aij;
            r_current = r_new;
            Metropolis=1;
            remove(i) = remove(i)+num_rewire;
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
    
    % if the target eigen-ratio is reached    
     if r_current <= min(r_target,r_optimal) %|| r_new <= min(r_target,r_optimal)
            OptimalNetwork=adj;
            r_target=r_current;%my code
            r_optimal=r_current;
    end
end

Summary.InitialLap=Laplacian2(AdjacencyTensor);
Summary.OptimalLap=Laplacian2(OptimalNetwork);
Summary.EigenRatio(NumberNode,1)=r_original;
Summary.EigenRatio(NumberNode,2)=r_optimal;
Summary.GeneralDegreeInitial{NumberNode}=CalculateDegree(AdjacencyTensor);
Summary.GeneralDegreeOptimal{NumberNode}=CalculateDegree(OptimalNetwork);
Summary.AdjTensorInitial{NumberNode}=AdjacencyTensor;
Summary.AdjTensorOptimal{NumberNode}=OptimalNetwork;
[Summary.AdjConvertedInitial{NumberNode},Summary.triangleInitial{NumberNode}]=ConvertTriangleToAdjacency(AdjacencyTensor);
[Summary.AdjConvertedOptimal{NumberNode},Summary.triangleOptimal{NumberNode}]=ConvertTriangleToAdjacency(OptimalNetwork);


%Plot the network
alpha=0.5;%transparence
shift=0.1;%shift=0;%0.018;%+(N-5)^1.2*(-0.0015);% avoid overlap
for Type=1:2%Initial and optimized network
    if Type==1
        triangle_index=Summary.triangleInitial{NumberNode};
        vtg=Summary.AdjConvertedInitial{NumberNode};
    elseif Type==2
        triangle_index=Summary.triangleOptimal{NumberNode};
        vtg=Summary.AdjConvertedOptimal{NumberNode};
    end
figure ('position', [00, 10, 500, 400]);
cc=jet(size(triangle_index,1));
for kk=1:size(triangle_index,1)
    AdjConverted=zeros(NumberNode,NumberNode);
    TempIndex=nchoosek(triangle_index(kk,:),2);
    for ll=1:size(TempIndex,1)
    AdjConverted(TempIndex(ll,1),TempIndex(ll,2)) = 1;
    AdjConverted(TempIndex(ll,2),TempIndex(ll,1)) = 1;%Symmetric
    end
drawCircGraphTriangle(AdjConverted,cc,kk,alpha,0,NumberNode); hold off;hold on;
%drawCircGraphTriangle2(AdjConverted,triangle_index(kk,:),cc,kk,alpha,shift,NumberNode); hold off;hold on;

end
figurenamehmm=[figureSubfolder,'\SecondOrder',num2str(NumberNode),'Type',num2str(Type),'.jpg'];
%saveas(gcf,figurenamehmm)
print(gcf,figurenamehmm,'-djpeg','-r300');

figure ('position', [00, 10, 500, 400]);
drawCircGraph(vtg)
figurenamehmm=[figureSubfolder,'\FirstOrder',num2str(NumberNode),'Type',num2str(Type),'.jpg'];
%saveas(gcf,figurenamehmm)
print(gcf,figurenamehmm,'-djpeg','-r300');
end
close all;

end

save(filename,'Summary');
end
%%
% function lap=Laplacian2(vtg)
% k_i = sum(sum(vtg,3),2)/2;
% k_ij = sum(vtg,3);
% lap = 2*diag(k_i) - k_ij;
% end


function lap=Laplacian1(vtg)
rowsum = sum(vtg,2);
lap = diag(rowsum) - vtg;
end

%% 
function GeneralDegree=CalculateDegree(OptimalNetwork)
for i=1:size(OptimalNetwork,1)
    GeneralDegree=(1/factorial(2))*sum(sum(OptimalNetwork,3),2);   
end
end

