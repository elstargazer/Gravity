ccc

G=6.67384e-11;
Resolution=0.25;
MaxDegreeTopo=300;
MaxDegreeGrav=70;
MaxTopoPower=10;
Rref=265000;
% rho=3457;
mu=17.28e9;
NTess=6;
aref=300000;
cref=300000;

refell=[aref cref];


%% SH shape model
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

[r1sh,lon,lat]=plm2xyz(lmcosi_shape,Resolution);

[lon,lat]=meshgrid(lon/180*pi,lat/180*pi);

[x1sh,y1sh,z1sh]=sph2cart(lon,lat,r1sh);

%% Tri shape model
TR=IcosahedronMesh;
TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;
xr=aref*TR_2.X(:,1);
yr=aref*TR_2.X(:,2);
zr=cref*TR_2.X(:,3);

s=size(FV);

[lambdar,fir,~]=cart2sph(xr,yr,zr);

r1tri=plm2xyz(lmcosi_shape,fir*180/pi,lambdar*180/pi);

[xr,yr,zr]=sph2cart(lambdar,fir,r1tri);

V=TetrahedronBobyVolume(xr',yr',zr',FV(:,1)',FV(:,2)',FV(:,3)');

% [vol,area] = triangulationVolume(TR_2,xr',yr',zr');
rho=mu/G/V;

%% Plotting shape models
figure; hold on; axis equal
tri_fig=trisurf(FV,xr,yr,zr,r1tri); shading interp
set(tri_fig,'edgecolor','k'); 
% trimesh_fig=trimesh(FV,x1,y1,z1,r1,'edgecolor','k','facecolor','none');

figure; hold on; axis equal
sh_fig=surf(x1sh,y1sh,z1sh,r1sh); shading interp


%% Reference surface

FiStep=1;
LambdaStep=1;
[lambdare,fire]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );

xre=aref*cos(lambdare).*cos(fire);
yre=aref*sin(lambdare).*cos(fire);
zre=cref*sin(fire);

%% gravity

lmcosi_gr=TopoSH2GravitySH(r1sh,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);
[gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_gr,xre,yre,zre);
[a_sh,~,~]=GravityComponents(gx,gy,gz,xre,yre,zre,aref,cref);

rho_tri=ones(1,s(1))*rho;

a_tri=RadGravityAccelerationTriDen(xr',yr',zr',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho_tri,refell);

%% Plotting residuals

res=(a_tri-a_sh)*1e5;

std(res(:))

% hist(res(:),20);

AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,res); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
% plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Residual anomaly','FontSize',20);
ylabel(cbar,'Acceletation [mGal]','FontSize',20);
% caxis(1e5*[min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])]);

