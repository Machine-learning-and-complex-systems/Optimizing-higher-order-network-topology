%A software package for optimizing synchronization of coupled oscillators with high-order networks
%(c) 2021 Ying Tang
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
InitialNet=1;%the two types of initialization as described in text: chose 1 or 2.
figurefolder=pwd;
figureSubfolder=[figurefolder,'\SecondOrder2'];
mkdir(figureSubfolder);
NumberNodes=[6, 10, 20,50,80,100];%Network sizes
%It requires the files generated from optiUndirected_rewire_high2.m and 
%need to specify the parameters matching that in optiUndirected_rewire_high2.m 
%to load the corresponding saved files. Please make sure the loaded files 
%are in the correct folder location.

for realization=1:10% the number of numerical replicates

filename=[figureSubfolder,'\Rewire2Scan4_',num2str(realization),'_InitialNet_',num2str(InitialNet),'.mat'];%Scan has target 1.8/1.4, Scan2 has target value 1
load(filename);

k=1;
for NumberNode=NumberNodes%, 50, 100]
lap2=Laplacian2(Summary.AdjTensorOptimal{NumberNode});
eigv2 = sort(real(eig(lap2)));
eigenratio.Second{k}(realization) = real(eigv2(end))/real(eigv2(2));

lap=Laplacian2(Summary.AdjConvertedOptimal{NumberNode});
eigv = sort(real(eig(lap)));
eigenratio.First{k}(realization) = real(eigv(end))/real(eigv(2));
k=k+1;
end
end


%% Eigenratios of second-order and corresponding first-order networks
for i=1:length(eigenratio.Second)
    Mean(i)=mean(eigenratio.Second{i});
    errlow(i)=std(eigenratio.Second{i});
    errhigh(i)=std(eigenratio.Second{i});
end

figure ('position', [00, 10, 500, 400]);
data = Mean';
bar(1:length(NumberNodes),data,'FaceColor',[0.5 0.5 0.5])                
hold on
er = errorbar(1:length(NumberNodes),data,errlow,errhigh,'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xticks(1:length(NumberNodes))
xticklabels(string(NumberNodes))
leg1 =ylabel('Eigenratios');
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
ylim([0, 10]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\EigenratioSecondOrder_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)

for i=1:length(eigenratio.First)
    Mean(i)=mean(eigenratio.First{i});
    errlow(i)=std(eigenratio.First{i});
    errhigh(i)=std(eigenratio.First{i});
end

figure ('position', [00, 10, 500, 400]);
data = Mean';
bar(1:length(NumberNodes),data,'FaceColor',[0.5 0.5 0.5])                
hold on
er = errorbar(1:length(NumberNodes),data,errlow,errhigh,'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xticks(1:length(NumberNodes))
xticklabels(string(NumberNodes))
leg1 =ylabel('Eigenratios');
ylim([0, 10]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
%ylim([0, 50]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\EigenratioFirstOrder_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)

