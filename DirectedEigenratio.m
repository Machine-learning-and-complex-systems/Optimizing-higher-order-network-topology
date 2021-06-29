clear;
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
lap=Laplacian2(AdjacencyTensor);
% find eigen-ratio of the original network
eigv = sort(eig(lap));
r_original = eigv(end)/eigv(2)

