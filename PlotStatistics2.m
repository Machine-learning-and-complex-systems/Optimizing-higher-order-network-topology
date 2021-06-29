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
figurefolder=pwd;
figureSubfolder=[figurefolder,'\ScanSecondOrder2'];
mkdir(figureSubfolder);
realization=1;
InitialNet=1;
filename=[figureSubfolder,'\Rewire2Scan4_',num2str(realization),'_InitialNet_',num2str(InitialNet),'.mat'];%Scan has target 1.8/1.4, Scan2 has target value 1
%filename='RewireScan2_2ndOrder_Indirect_InitialNet_4';% First need to generate the file from optiUndirected-rewire-high.m
%filename2=[pwd,'\',filename,'.mat'];
load(filename);
Ntotal=100;%Network size


%% Degree distribution
figure ('position', [00, 10, 500, 400]);
Max=max(max(Summary.EigenRatio));
cc=parula(Ntotal);
plot([0:Max+1],[0:Max+1],'k--','linewidth',1.5);hold on;
for n=1:Ntotal
    if Summary.EigenRatio(n,1)~=0 && Summary.EigenRatio(n,2)~=0
    plot(Summary.EigenRatio(n,1),Summary.EigenRatio(n,2),'o','MarkerSize',10,...
    'MarkerEdgeColor',cc(n,:),'MarkerFaceColor',cc(n,:));
    hold on;
    end
end
legend('The line of equal value','location','north');
%set(gca, 'YScale', 'log');ylim([1e-3 1e0]);
colormap(parula(Ntotal))
hc = colorbar;
leg1=ylabel(hc, 'Network size $N$')
set(leg1,'Interpreter','latex');
cb = [1,linspace(10,Ntotal,10)];
cb1 = linspace(0,1,11);
set(hc, 'YTick',cb1, 'YTickLabel',cb)
xlim([0 Max]);ylim([0 Max]);
 xlim([0 15]);ylim([0 15]);
leg1=ylabel('Eigenratio of rewired network');
set(leg1,'Interpreter','latex');
leg1=xlabel('Eigenratio of initial network');
set(leg1,'Interpreter','latex');
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\Eigenratio_size',num2str(Ntotal),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)


%% Degree distribution
%It requires the files generated from optiUndirected_rewire_high2.m and 
%need to specify the parameters matching that in optiUndirected_rewire_high2.m 
%to load the corresponding saved files. Please make sure the loaded files 
%are in the correct folder location.
for nn=[6,10:5:100]%Network sizes
   % h1 =histogram(Summary.GeneralDegreeInitial{nn});
   % h2 =histogram(Summary.GeneralDegreeOptimal{nn});
    
figure ('position', [00, 10, 500, 400]);
%edges = [0:5:100];
width=1;
%
%hold on;

%h1.Normalization = 'probability';
%h1.BinWidth = width;
%h2.Normalization = 'probability';
%h2.BinWidth = width;
y = [Summary.GeneralDegreeInitial{nn},...
    Summary.GeneralDegreeOptimal{nn}];
h=hist(y,0:12);
bar(h/max(sum(h,1)));
legend('initial','rewired','location','northeast');
%set(gca, 'YScale', 'log');ylim([1e-3 1e0]);
xlim([1 12]);
leg1 =ylabel('Distribution $P(k_i^{(2)})$');
set(leg1,'Interpreter','latex');
leg2 = xlabel('Generalized degree $k_i^{(2)}$');
set(leg2,'Interpreter','latex');
%ylim([0, 1]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\DegreeDist_size',num2str(nn),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)
close all;
end
