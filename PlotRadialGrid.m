function PlotRadialGrid(ri)

s=size(ri);

fi=linspace(90,-90,s(1))/180*pi;
lambda=linspace(0,360,s(2))/180*pi;

[lambdai,fii]=meshgrid(lambda,fi);

[xi,yi,zi]=sph2cart(lambdai,fii,ri);

surf(xi,yi,zi,ri);

 lighting phong
 shading interp
 axis tight equal off
 light('Position',[-1 1 0],'Style','infinite');

% [fig,ax]=AGUaxes;
% 
% surfm(fii,lambdai,ri);
% 
% colorbar