% Draw a circular representation of a graph with the nodes ordered by degree
% Strategy: position vertices in a regular n-polygon
%
% INPUTs: adj, nxn - adjacency matrix
% OUTPUTs: plot
%
% Other routines used: degrees.m
% GB: last updated, Nov 29 2012

function [] = drawCircGraphTriangle(adj,adj2,cc,kk,alpha,shift,N)

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
  plot1=plot(x(k)*(1+kk*shift),y(k)*(1+kk*shift),'ko','markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');%'markersize',30)
  text(1.15*x(k),1.15*y(k),strcat(num2str(2*k-1)),"fontsize",20);
  plot1.Color(4) = alpha;
  hold off; hold on;
end

edges=find(adj>0);
for e=1:length(edges)
    [ii,jj]=ind2sub([n,n],edges(e));
    %disp(e);
    nodeX(e)=x(ii);
    nodeY(e)=y(ii);
    %lh.Color = [lh.Color alpha];
    hold off; hold on;
end


edges=find(adj2>0);
for e=1:length(edges)
    [ii,jj]=ind2sub([n,n],edges(e));
    plot_arrow( x(ii),y(ii),(x(ii)+x(jj))/2,(y(ii)+y(jj))/2,'linewidth',1,'color',cc(kk,:),'facecolor',cc(kk,:), 'edgecolor',cc(kk,:),'headwidth',0.2,'headheight',0.2);
end


%nodeX=unique(nodeX);
%nodeY=unique(nodeY);
%assignin('base','nodeY',nodeY);
%assignin('base','nodeX',nodeX);
%plot(intersect(pgon{1},pgon{2}),'EdgeColor','none','FaceColor','m')
%Plot shaded triangles:
w=10; %width of triangle
ar=0.866; % Aspect ratio for equilateral triangle
h=ar*w;%height of triangle
p1=patch(nodeX,nodeY,cc(kk,:), 'facecolor',cc(kk,:), 'facealpha',alpha,'linewidth',1, 'edgecolor',cc(kk,:));%, 'edgealpha',alpha);%plotting triangle in white color
%p1.FaceAlpha = alpha;%'interp' ; 
%p1.FaceVertexAlphaData = alpha;    % Set constant transparency 
%p1.FaceAlpha = 'flat' ;  
set(gca,'xticklabel',[],'yticklabel',[])%setting background black and removing labels
%daspect([1 1 1]);%equal data unit length along x and y axis


edg=1.15;%1+N*0.1;
axis([-edg edg -edg edg])
axis square
axis off
