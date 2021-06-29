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
figureSubfolder=[figurefolder,'\DirectedNetworkAsymmetry2'];
mkdir(figureSubfolder);
NumberNodes=[6, 10, 20,50,80,100];%Network sizes
p=0.25; %percentile for the bar graph errorbar

for realization=1:30%0% the number of numerical replicates

filename2=[figureSubfolder,'\RemoveScan_',num2str(realization),'_InitialNet_',num2str(InitialNet),'.mat'];%Scan has target 1.8/1.4, Scan2 has target value 1
load(filename2);

k=1;
for NumberNode=NumberNodes%, 50, 100]
lap2=Laplacian2(Summary.AdjTensorOptimal{NumberNode});
eigv2 = sort(real(eig(lap2)));
eigenratio.Second{k}(realization) = real(eigv2(end))/real(eigv2(2));
AsymmetryIndex1{k}(realization) = mean(Summary.AsymmetryIndex1{NumberNode});
AsymmetryIndex2{k}(realization) = mean(Summary.AsymmetryIndex2{NumberNode});

k=k+1;
end
end


%% Violin plots


colors = hsv(length(AsymmetryIndex1));
figure ('position', [00, 10, 500, 400]);
for ii=1:length(AsymmetryIndex1)
violinPlot(AsymmetryIndex1{ii}','xValues',ii, 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 0,...
    'color',  mat2cell(colors(ii, : ), 1)); hold on;
violinPlot(AsymmetryIndex1{ii}','xValues',ii, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0,...
    'color',  mat2cell(colors(ii, : ), 1));
end
alpha(.5);
set(gca, 'xlim', [0 length(AsymmetryIndex1)+1],'XTick',1:length(AsymmetryIndex1),'XTickLabel',string(NumberNodes));
set(gca,'TickLabelInterpreter','none');
box on;
%xtickangle(45)
leg1 =ylabel({'Average asymmetric measure\\','for each node'});
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
ylim([0, 10]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\ViolinAsymmetry1Directed_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)


colors = hsv(length(AsymmetryIndex2));
figure ('position', [00, 10, 500, 400]);
for ii=1:length(AsymmetryIndex2)
violinPlot(AsymmetryIndex2{ii}','xValues',ii, 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 0,...
    'color',  mat2cell(colors(ii, : ), 1)); hold on;
violinPlot(AsymmetryIndex2{ii}','xValues',ii, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0,...
    'color',  mat2cell(colors(ii, : ), 1));
end
alpha(.5);
set(gca, 'xlim', [0 length(AsymmetryIndex2)+1],'XTick',1:length(AsymmetryIndex2),'XTickLabel',string(NumberNodes));
set(gca,'TickLabelInterpreter','none');
box on;
%xtickangle(45)
leg1 =ylabel({'Average asymmetric measure\\','for each node'});
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
ylim([0, 25]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\ViolinAsymmetry2Directed_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)

%% Asymmetry1 of second-order directed networks
for i=1:length(eigenratio.Second)
    Mean(i)=mean(AsymmetryIndex1{i});
%     errlow(i)=std(AsymmetryIndex1{i});
%     errhigh(i)=std(AsymmetryIndex1{i});
    errlow(i)=prctile(AsymmetryIndex1{i},p); 
     errhigh(i)=prctile(AsymmetryIndex1{i},p); 
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
leg1 =ylabel({'Average asymmetric measure\\','for each node'});
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
ylim([0, 10]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\Asymmetry1Directed_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)




%% Asymmetry2 of second-order directed networks
for i=1:length(eigenratio.Second)
    Mean(i)=mean(AsymmetryIndex2{i});
%     errlow(i)=std(AsymmetryIndex2{i});
%     errhigh(i)=std(AsymmetryIndex2{i});
     errlow(i)=prctile(AsymmetryIndex1{i},p); 
     errhigh(i)=prctile(AsymmetryIndex1{i},p); 
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
leg1 =ylabel({'Average asymmetric measure\\','for each node'});
set(leg1,'Interpreter','latex');
leg2 = xlabel('Network size');
set(leg2,'Interpreter','latex');
ylim([0, 20]);
%set(gca, 'YScale', 'log');ylim([0.9 1e3]);
set(gca,'FontSize',16,'linewidth',1.5);
figurenamehmm=[figureSubfolder,'\Asymmetry2Directed_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)
%% Eigenratios of second-order directed networks
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
figurenamehmm=[figureSubfolder,'\EigenratioDirected_InitialNet_',num2str(InitialNet),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm)
