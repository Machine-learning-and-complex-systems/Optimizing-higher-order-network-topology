% Draw a circular representation of a graph with the nodes ordered by degree
% Strategy: position vertices in a regular n-polygon
%
% INPUTs: adj, nxn - adjacency matrix
% OUTPUTs: plot
%
% Other routines used: degrees.m
% GB: last updated, Nov 29 2012

function [] = drawCircGraph(adj)

n = size(adj,1); % number of nodes
[~, Y] = sort(degrees(adj));   % Y - sorted nodal indices
angl = 2*pi/n; % rotation angle

for k=1:n
  x(Y(k)) = real(exp(angl*(k-1)*i));
  y(Y(k)) = imag(exp(angl*(k-1)*i));
end

% regulate the font size with respect to graph size
if n<=50;
  fontsize = 10;
elseif n>50 && n<=150
  fontsize = 8;
else
  fontsize = 6;
end

for k=1:n
  plot(x(k),y(k),'ko','markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');%'markersize',30)
  %text(1.15*x(k),1.15*y(k),strcat('v',num2str(k)),"fontsize",fontsize);
  hold off; hold on;
end

edges=find(adj>0);
set(gcf,'Color',[1,1,1])

c=jet(6);
for e=1:length(edges)
    [ii,jj]=ind2sub([n,n],edges(e));
    %line([x(ii) x(jj)],[y(ii) y(jj)],'Color',[0.000000 0.50000 1.000000], 'LineWidth',2);
    line([x(ii) x(jj)],[y(ii) y(jj)],'Color','b', 'LineWidth',1);
    hold off; hold on;
end
edg=1.15;%1+N*0.1;
axis([-edg edg -edg edg])
axis square
axis off
