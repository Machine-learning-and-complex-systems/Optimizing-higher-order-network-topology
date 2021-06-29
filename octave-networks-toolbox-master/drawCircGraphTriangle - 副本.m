% Draw a circular representation of a graph with the nodes ordered by degree
% Strategy: position vertices in a regular n-polygon
%
% INPUTs: adj, nxn - adjacency matrix
% OUTPUTs: plot
%
% Other routines used: degrees.m
% GB: last updated, Nov 29 2012

function [] = drawCircGraphTriangle(adj,cc,kk,alpha,shift,N)

n = size(adj,1); % number of nodes
[~, Y] = sort(degrees(adj));   % Y - sorted nodal indices
%disp(Y);
angl = 2*pi/n; % rotation angle

for k=1:n
%   x(Y(k)) = real(exp(angl*(k-1)*i));
%   y(Y(k)) = imag(exp(angl*(k-1)*i));
  x(k) = real(exp(angl*(k-1)*i));
  y(k) = imag(exp(angl*(k-1)*i));
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
  plot(x(k)*(1+kk*shift),y(k)*(1+kk*shift),'ko','markersize',15,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);%'markersize',30)
  %text(1.15*x(k),1.15*y(k),strcat('v',num2str(k)),"fontsize",fontsize);
  hold off; hold on;
end

edges=find(adj>0);
set(gcf,'Color',[1,1,1])

for e=1:length(edges)
    [ii,jj]=ind2sub([n,n],edges(e));
    lh=line([x(ii)*(1+kk*shift) x(jj)*(1+kk*shift) ],[y(ii)*(1+kk*shift) y(jj)*(1+kk*shift)],'Color',cc(kk,:), 'LineWidth',2);
    lh.Color = [lh.Color alpha];
    hold off; hold on;
end
edg=1+N*0.1;
axis([-edg edg -edg edg])
axis square
axis off
