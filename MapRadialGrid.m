function [fig,ax]=MapRadialGrid(varargin)


ri=varargin{1};

if (nargin==1)
s=size(ri);

fi=linspace(-90,90,s(1))/180*pi;
lambda=linspace(0,360,s(2))/180*pi;

[lambdai,fii]=meshgrid(lambda,fi);

% [xi,yi,zi]=sph2cart(lambdai,fii,ri);

% surf(xi,yi,zi,ri);

%  lighting phong
%  shading interp
%  axis tight equal off
%  light('Position',[-1 1 0],'Style','infinite');

[fig,ax]=AGUaxes;


set(ax,'fontsize',20);

minr=min(min(ri))
maxr=max(max(ri))

StepSize=5;

surfm(fii,lambdai,ri);
% StandardLight
colorbar

% Labels=floor(minr/StepSize)*StepSize:StepSize:floor(maxr/StepSize)*StepSize+StepSize;
% 
% [c,h]=contourm(fii,lambdai,ri,Labels,'LineWidth',3);
% 
 xlabel(['min = ' num2str(minr) '; max = ' num2str(maxr)],'fontsize',25);
 
 
% 
% ht = clabelm(c,h);
% 
% set(ht,'Color','k','BackgroundColor','white','FontWeight','bold','FontSize',20)

% colorbar('peer',ax);

else 
    
ri=varargin{1};
minr=varargin{2};   
maxr=varargin{3};   

s=size(ri);

fi=linspace(-90,90,s(1))/180*pi;
lambda=linspace(0,360,s(2))/180*pi;

[lambdai,fii]=meshgrid(lambda,fi);

% [xi,yi,zi]=sph2cart(lambdai,fii,ri);

% surf(xi,yi,zi,ri);

%  lighting phong
%  shading interp
%  axis tight equal off
%  light('Position',[-1 1 0],'Style','infinite');

[fig,ax]=AGUaxes;

caxis([minr maxr])

set(ax,'fontsize',20);

StepSize=5;

surfm(fii,lambdai,ri);
% StandardLight
colorbar

% Labels=floor(minr/StepSize)*StepSize:StepSize:floor(maxr/StepSize)*StepSize+StepSize;
% 
% [c,h]=contourm(fii,lambdai,ri,Labels,'LineWidth',3);
% 
 xlabel(['min = ' num2str(minr) '; max = ' num2str(maxr)],'fontsize',25);
 
 
% 
% ht = clabelm(c,h);
% 
% set(ht,'Color','k','BackgroundColor','white','FontWeight','bold','FontSize',20)

% colorbar('peer',ax);   
    
    
end