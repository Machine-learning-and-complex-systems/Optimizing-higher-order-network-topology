function [AdjConverted,triangle_index]=ConvertTriangleToAdjacency(AdjacencyTensor)


AdjConverted=zeros(size(AdjacencyTensor,1),size(AdjacencyTensor,1));
idx = find(AdjacencyTensor == 1);
[i1, i2, i3] = ind2sub(size(AdjacencyTensor), idx);
triangle_index=[i1(1:length(idx)), i2(1:length(idx)), i3(1:length(idx))];
triangle_index = unique(sort(triangle_index,2),'rows');

for kk=1:size(triangle_index,1)
    %triangle_links=perms(triangle_index(kk,:));
    TempIndex=nchoosek(triangle_index(kk,:),2);
    for ll=1:size(TempIndex,1)
    AdjConverted(TempIndex(ll,1),TempIndex(ll,2)) = 1;
    AdjConverted(TempIndex(ll,2),TempIndex(ll,1)) = 1;%Symmetric
    end
end
for kk=1:size(AdjConverted,1)
    AdjConverted(kk,kk)=0;
end

end