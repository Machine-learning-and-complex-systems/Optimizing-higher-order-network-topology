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
clc;
addpath([pwd,'\octave-networks-toolbox-master']);
%% Directed network counterexample
% NumberNode=7;
% for j=1:NumberNode%NumberNode
%     id=mod([1,3,4]+(j-1),NumberNode)+1;
%     P = perms(id); %symmetric tensor for each triangle
%     for k=1:size(P,1)
%     AdjacencyTensor(P(k,1),P(k,2),P(k,3))=1;
%     end
% end
% %Triangles={[1,2,3],[1,2,5],[1,3,4],[1,4,5],[2,3,6],[2,5,6],[3,4,6],[4,5,6]};% initialize 8 triangles for Shi
% %Triangles={[1,2,3],[1,2,6],[1,3,5],[1,4,6],[2,3,4],[2,4,5],[3,5,6],[4,5,6]};% initialize 8 triangles for Shi
% % for j=1:length(Triangles)
% %     id=Triangles{j};
% % P = perms(id); %symmetric tensor for each triangle
% % for k=1:size(P,1)
% % AdjacencyTensor(P(k,1),P(k,2),P(k,3))=1;
% % end
% % end

NumberNode=7;
for j=1:NumberNode%NumberNode
    id=mod([1,2,3]+(j-1),NumberNode)+1;
    P = perms(id); %symmetric tensor for each triangle
    for k=1:size(P,1)
    AdjacencyTensor(P(k,1),P(k,2),P(k,3))=0;
    end
end

DirectedTriangles=...
   {[1,2,4], [2,3,5], [3,4,6], [4,5,7], [5,6,1], [6,7,2], [7,1,3],...
    [1,3,5], [2,4,6], [3,5,7], [4,6,1], [5,7,2], [6,1,3], [7,2,4],...
    [1,6,7], [2,7,1], [3,1,2], [4,2,3], [5,3,4], [6,4,5], [7,5,6]};%,...
    %[1,2,5], [2,3,6], [3,4,7], [4,5,1], [5,6,2], [6,7,3], [7,1,4]};%,...
   % [1,2,6], [2,3,7], [3,4,1], [4,5,2], [5,6,3], [6,7,4], [7,1,5]}; 
for i=1:length(DirectedTriangles)
    AdjacencyTensor(DirectedTriangles{i}(1),DirectedTriangles{i}(2),DirectedTriangles{i}(3))=1;
    AdjacencyTensor(DirectedTriangles{i}(1),DirectedTriangles{i}(3),DirectedTriangles{i}(2))=1;
   %AdjacencyTensor(DirectedTriangles{i}(2),DirectedTriangles{i}(1),DirectedTriangles{i}(3))=1;
    %AdjacencyTensor(DirectedTriangles{i}(2),DirectedTriangles{i}(3),DirectedTriangles{i}(1))=1;
   %AdjacencyTensor(DirectedTriangles{i}(3),DirectedTriangles{i}(1),DirectedTriangles{i}(2))=1;
    %AdjacencyTensor(DirectedTriangles{i}(3),DirectedTriangles{i}(2),DirectedTriangles{i}(1))=1;
end

[Summary.AdjConvertedInitial{NumberNode},Summary.triangleInitial{NumberNode}]=ConvertTriangleToAdjacency(AdjacencyTensor);

for N=NumberNode%5:10
alpha=0.5;%transparence
%shift=0;%0.018;%+(N-5)^1.2*(-0.0015);% avoid overlap
shift=0.1;%+(N-5)^1.2*(-0.0015);% avoid overlap
clear triangle_index;

triangle_index=Summary.triangleInitial{N};
figure ('position', [00, 10, 500, 400]);
cc=jet(size(triangle_index,1));
for kk=1:size(triangle_index,1)
    AdjConverted=zeros(N,N);
    TempIndex=nchoosek(triangle_index(kk,:),2);
    for ll=1:size(TempIndex,1)
    AdjConverted(TempIndex(ll,1),TempIndex(ll,2)) = 1;
    AdjConverted(TempIndex(ll,2),TempIndex(ll,1)) = 1;%Symmetric
    end
    %disp(AdjConverted);
%drawCircGraphTriangle(AdjConverted,cc,kk,alpha,shift,N); hold off;hold on;
drawCircGraphTriangle2(AdjConverted,triangle_index(kk,:),cc,kk,alpha,shift,N); hold off;hold on;

end
figurenamehmm=[pwd,'\New2SecondOrder',num2str(N),'.jpg'];
print(gcf,figurenamehmm,'-djpeg','-r300');


end



%% Figure 1
% vtg=[0,1,1,1,1;...
%      1,0,1,1,1;...
%      1,1,0,1,1;...
%      1,1,1,0,1;...
%      1,1,1,1,0];
vtg=[0,1,0;...
     0,0,1;...
     1,0,0];

%  vtg=[0,1;...
%      0,0];

%   vtg=[1];
figure ('position', [00, 10, 500, 400]);
drawCircGraph(vtg)
figurenamehmm=[pwd,'\SchematicNetwork_order',num2str(size(vtg,1)),'.jpg'];
%saveas(gcf,figurenamehmm)
print(gcf,figurenamehmm,'-djpeg','-r300');

% figure ('position', [00, 10, 500, 400]);
% vtg=Summary.AdjConvertedOptimal{N};
% drawCircGraph(vtg)
% figurenamehmm=[figureSubfolder,'\Network1stOrderOptimal_',filename,'_size',num2str(N),'.jpg'];
% saveas(gcf,figurenamehmm)