% ccc
close all
%% Input parameters

rho_mean=3457;
Rref=265101;

xcm=0;
ycm=0;
zcm=0;

%% Make tri shape model
NTess=1;
% filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
% [FV,xi,yi,zi]=Grid2TriShapeModel(filename,NTess);

Radius=265000;
TR=IcosahedronMesh;

TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;

xi=Radius*TR_2.X(:,1);
yi=Radius*TR_2.X(:,2);
zi=Radius*TR_2.X(:,3);

xi=xi+xcm;
yi=yi+xcm;
zi=zi+zcm;

[lambdai,fii,ri]=cart2sph(xi,yi,zi);

ri=sqrt(xi.*xi+yi.*yi+zi.*zi);

% density model

s=size(FV);
% rho=rho_mean+0.01*rho_mean*randn(1,s(1));
rho=ones(1,s(1))*rho_mean;

N=55;

rho(N)=rho(N)*8;

xc=(xi(FV(N,1))+xi(FV(N,2))+xi(FV(N,3)))/3;
yc=(yi(FV(N,1))+yi(FV(N,2))+yi(FV(N,3)))/3;
zc=(zi(FV(N,1))+zi(FV(N,2))+zi(FV(N,3)))/3;

[lambdac,fic,~]=cart2sph(xc,yc,zc);



%% Reference surface

FiStep=5;
LambdaStep=5;

minfi=fic-25/180*pi
maxfi=fic+25/180*pi

minlambda=lambdac-25/180*pi
maxlambda=lambdac+25/180*pi



[lambdar,fir]=meshgrid(minlambda:LambdaStep/180*pi:maxlambda,...
    minfi:FiStep/180*pi:maxfi);


xr=Rref*cos(lambdar).*cos(fir);
yr=Rref*sin(lambdar).*cos(fir);
zr=Rref*sin(fir);

rr=sqrt(xr.*xr+yr.*yr+zr.*zr);

%% Computing acceleration


Ver1=FV(:,1);
Ver2=FV(:,2);
Ver3=FV(:,3);


% Computing Potential and Accelertion
% [U,ax,ay,az]=GravityPotentialAccelerationTri(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,ro);

[ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
[U_m,ax_m,ay_m,az_m]=GravityPotentialAccelerationTri(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_mean);
U=GravityPotentialTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);

% Computing Accelertion derivatives

% d=AccelerationDerivativeTri(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,ro);


% a=sqrt(ax.*ax+ay.*ay+az.*az);
% a_m=sqrt(ax_m.*ax_m+ay_m.*ay_m+az_m.*az_m);
% a_diff=sqrt((ax-ax_m).^2+(ay-ay_m).^2+(az-az_m).^2)
% % [min(a_diff(:)) max(a_diff(:))]
% [min(a(:)) max(a(:))
% [min(a_m(:)) max(a_m(:))]
 

[fig,axes1]=AGUaxes;
pcolorm(fir,lambdar,U); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
colorbar
% set(axes1,'MapLatLimit',[minfi maxfi],'MapLatLimit',[minlambda maxlambda]);


plotm(fic,lambdac,'bo','MarkerFaceColor','k','MarkerSize',9);

%% Plotting results

figure('Position',[1 1 1200 1200]); hold on; axis equal
tri_fig=trisurf(FV,xi,yi,zi,ri); shading interp
set(tri_fig,'edgecolor','k','facecolor','none'); 


% trimesh_fig=trimesh(FV,xi,yi,zi,ri,'edgecolor','k');

% surf(xr,yr,zr,a); 

 plot3(xr,yr,zr,'.k','MarkerSize',1);
 quiver3(xr,yr,zr,-(ax-ax_m),-(ay-ay_m),-(az-az_m),5,'MarkerSize',1,'LineWidth',2)
 axis(1.1*[-300000 300000 -300000 300000 -300000 300000])
 