ccc

%% Input values

% input gravity and shape model files
GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';

NTess=5;
NTessPlus=0;
k=1;

if (k<-NTessPlus)
    error('Not enough data points');
end


N_trunc=15;

% mantle
a2=257000;
c2=207000;

% core
a3=117000;
c3=105000;

% reference ellipsoid
aref=295000;
cref=295000;

% mean density
rho_mean=3457;

% internal structure densities
rho_core=7400;
rho_mantle=3160;
rho_crust=2885;

% core offset
% x3c=2.2291e+03;
% y3c=9.4132e+03;
% z3c=-29.7072;

%% Load gravity model

[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);
lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);


lmcosi_pot=plm2pot(AddZeroHarm(lmcosi_obs,1),aref,mu,Rref,1,'nothing');

%% Load topography model

[FV,x1,y1,z1]=Grid2TriShapeModel(filename,NTess);
x1=x1*1000;
y1=y1*1000;
z1=z1*1000;

s=size(FV);
[lambdai,fii,r1]=cart2sph(x1,y1,z1);

% MaxDegreeTopo=7;
% load VestaHASTALAVESTAshape_sh720.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% 
% r1sh=plm2xyz(lmcosi_shape,fii*180/pi,lambdai*180/pi);
% [x1,y1,z1]=sph2cart(lambdai,fii,r1sh);

x1m=(x1(FV(:,1))+x1(FV(:,2))+x1(FV(:,3)))/3;
y1m=(y1(FV(:,1))+y1(FV(:,2))+y1(FV(:,3)))/3;
z1m=(z1(FV(:,1))+z1(FV(:,2))+z1(FV(:,3)))/3;

%% Create reference surface

refell=[aref aref];

TR=IcosahedronMesh;
TR_2r=SubdivideSphericalMesh(TR,NTess+NTessPlus); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FVr=TR_2r.Triangulation;
xr=aref*TR_2r.X(:,1);
yr=aref*TR_2r.X(:,2);
zr=cref*TR_2r.X(:,3);
xr=(xr(FVr(:,1))+xr(FVr(:,2))+xr(FVr(:,3)))/3;
yr=(yr(FVr(:,1))+yr(FVr(:,2))+yr(FVr(:,3)))/3;
zr=(zr(FVr(:,1))+zr(FVr(:,2))+zr(FVr(:,3)))/3;

% amp=5000;
% xr=xr+amp*rand(size(xr));
% yr=yr+amp*rand(size(yr));
% zr=zr+amp*rand(size(zr));

[lambdar,fir,~]=cart2sph(xr,yr,zr);
[xr,yr,zr]=sph2cart(lambdar,fir,1);

xr=aref*xr;
yr=aref*yr;
zr=cref*zr;

[~,~,rr]=cart2sph(xr,yr,zr);



%% Plot

% figure; hold on; axis equal
% tri_fig=trisurf(FV,x1,y1,z1,r1); shading interp
% set(tri_fig,'edgecolor','k'); 
% trimesh_fig=trimesh(FV,x1,y1,z1,r1,'edgecolor','k','facecolor','none');
% 
% % figure; hold on; axis equal
% % tri_fig=trisurf(FV,x1sh,y1sh,z1sh,r1sh); shading interp
% % set(tri_fig,'edgecolor','k'); 
% % trimesh_fig=trimesh(FV,x1sh,y1sh,z1sh,r1sh,'edgecolor','k','facecolor','none');
% 
% % surf(xr,yr,zr,a); 
% 
% plot3(xr,yr,zr,'.k','MarkerSize',10);
% axis(1.4*[-300000 300000 -300000 300000 -300000 300000])
% alpha(0.2);

%% Constructing internal structure

x2=a2*cos(lambdai).*cos(fii);
y2=a2*sin(lambdai).*cos(fii);
z2=c2*sin(fii);

x3=a3*cos(lambdai).*cos(fii);
y3=a3*sin(lambdai).*cos(fii);
z3=c3*sin(fii);

% figure; hold on; axis equal
% tri_fig1=trisurf(FV,x1,y1,z1,r1); 
% tri_fig2=trisurf(FV,x2,y2,z2); 
% tri_fig3=trisurf(FV,x3,y3,z3); 
% 
% shading interp; alpha(0.5);

rho1=ones(1,s(1))*rho_crust;
rho2=ones(1,s(1))*rho_mantle-rho1;
rho3=ones(1,s(1))*rho_core-rho2-rho1;


%% Computing gravity
% [gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_obs,xr,yr,zr);
% [g_up,g_east,g_north]=GravityComponents(gx,gy,gz,xr,yr,zr,aref,cref);

% U_o=GravityPotential(mu,Rref,lmcosi_obs,xr,yr,zr);

% generating observed potential
% rho1_i=ones(1,s(1))*rho_crust; 
% rho1_i=rho1_i+200*randn(size(rho1_i));
% rho2_i=ones(1,s(1))*rho_mantle-rho1_i;
% rho3_i=ones(1,s(1))*rho_core-rho2_i-rho1_i;
% 
% 
% U1_o=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho1_i);
% U2_o=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho2_i);
% U3_o=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho3_i);
% 
% U_o=U1_o+U2_o+U3_o;

tic
[U_o,~,~,~]=plm2xyz(lmcosi_pot,fir*180/pi,lambdar*180/pi);
t_pot_obs=toc

% gtot=sqrt(gx.*gx+gy.*gy+gz.*gz);

% [ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% atot=sqrt(ax.*ax+ay.*ay+az.*az);
% U_c=GravityPotentialTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% d2=DeltaSqGravity(rho,xi,yi,zi,FV,xr,yr,zr,gx,gy,gz)

%% Numerical Jacobian

% Jacobian linear
% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho,refell);
% a_up0=a_up;

tic
% U1_c=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho1);
% U2_c=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho2);
% U3_c=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho3);


% Potential from triaxial ellipsoid

t_pot_calc=toc

U_c=U1_c+U2_c+U3_c;

U_c0=U_c;

delta_rho=50;
Nelems=numel(rho1);
Nobs=numel(xr);


rho1=ones(1,s(1))*rho_crust;
rho2=ones(1,s(1))*rho_mantle-rho1;
rho3=ones(1,s(1))*rho_core-rho2-rho1;

% Jacobian numerical
tic;
J1=GravityPotentialTriDenDer(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);
J2=GravityPotentialTriDenDer(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);
t_deriv=toc;

jac=J1-J2;

%% Inversion

L=U_o-U_c; % initial residuals

sumsq=sum(L.^2);

jac_orig=jac;
% jac=jac_orig;


% reducing jacobian

kp=4^k;

if (k>0)    
    jac_red=zeros(Nobs,Nelems/kp);
    for n=1:Nelems/kp
        jac_red(:,n)=sum(jac(:,n:Nelems/kp:Nelems),2);
    end    
    jac=jac_red;
end


A=(jac'*jac);
b=jac'*L;

drho1=A\b;

std(drho1)



% drho1=((jac'*jac)^-1)*jac'*L;

 range_crust=500;
 
% Constrained inversion

% lb=zeros(size(rho1))-range_crust;
% ub=zeros(size(rho1))+range_crust;
% drho1 = lsqlin(A,b,[],[],[],[],lb,ub);

% rho1_s=rho1+drho1';


for n=1:Nelems/kp
rho1_s(n:Nelems/kp:Nelems)=rho1(n:Nelems/kp:Nelems)+drho1(n);
end




rho2_s=ones(1,s(1))*rho_mantle-rho1_s;
rho3_s=ones(1,s(1))*rho_core-rho2_s-rho1_s;

% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho,refell);
U1_c=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho1_s);
U2_c=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho2_s);
U3_c=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho3_s);

U_c=U1_c+U2_c+U3_c;

figure; hold on;
set(gca,'FontSize',20);
% plot((a_up-g_up)./g_up*100)
plot((L)./U_o*100,'.b');
plot((U_c-U_o)./U_o*100,'.r');
legend({'Initial residuals','Final residuals'},'FontSize',20);
ylabel('Potential residual [%] ','FontSize',20);



figure; hold on; 
set(gca,'FontSize',20);
plot(rho1_s,'.');
plot(ones(size(rho1_s))*rho_crust+range_crust,'--k');
plot(ones(size(rho1_s))*rho_crust-range_crust,'--k');
plot(ones(size(rho1_s))*rho_crust,'-k');
ylabel('Density [kg/m^3]','FontSize',20);


%% Plot Resuduals

FiStep=2;
LambdaStep=2;
[lambdare,fire]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );

xre=aref*cos(lambdare).*cos(fire);
yre=aref*sin(lambdare).*cos(fire);
zre=cref*sin(fire);

rre=sqrt(xre.*xre+yre.*yre+zre.*zre);

% Plotting final residuals

% a_up_re=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho,refell);
U1_c_re=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho1_s);
U2_c_re=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho2_s);
U3_c_re=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho3_s);

U_c_re=U1_c_re+U2_c_re+U3_c_re;

% [gx_re,gy_re,gz_re]=GravityAcceleration(mu,Rref,lmcosi_obs,xre,yre,zre);
% [g_up_re,g_east_re,g_north_re]=GravityComponents(gx_re,gy_re,gz_re,xre,yre,zre,aref,cref);

% U_o_re=GravityPotential(mu,Rref,lmcosi_obs,xre,yre,zre);

[U_o_re,~,~,~]=plm2xyz(lmcosi_pot,fire(:)*180/pi,lambdare(:)*180/pi);
U_o_re=reshape(U_o_re,size(fire));



%% Comparing initial potential (for fake generated initial potential)
% U1_o_re=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho1_i);
% U2_o_re=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho2_i);
% U3_o_re=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho3_i);
% 
% U_o_re=U1_o_re+U2_o_re+U3_o_re;


U1_init=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho1);
U2_init=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho2);
U3_init=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho3);

U_c_init=U1_init+U2_init+U3_init;


U_res_init=U_o_re-U_c_init;
 

%% Initial Potential residuals
% [U_o_re2,lon,lat,Plm]=plm2xyz(lmcosi_pot,fire*180/pi,lambdare*180/pi);
AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,U_res_init); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
% plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Initial residuals','FontSize',20);
ylabel(cbar,'Potential [m^{2} s^{-2} ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
% caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

%% Plotting observed potential 


AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,U_o_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from SH','FontSize',20);
ylabel(cbar,'Potential [m^2 s^-2 ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])


%% Plotting calculated potential
AGUaxes
% pcolorm(fire,lambdare,a_up_re); shading interp;
pcolorm(fire,lambdare,U_c_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from polyhedra','FontSize',20);
ylabel(cbar,'Potential [m^2 s^-2 ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

% res=(a_up_re-g_up_re)*1e5;
U_res_final=(U_c_re-U_o_re);

figure; 
set(gca,'FontSize',20);
hist(U_res_final(:)./U_o_re(:)*100,20);
xlabel('Residual potential [%]','FontSize',20);

%% Final potential resuduals
AGUaxes
pcolorm(fire,lambdare,U_res_final); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Residual potential [m^2 s^-2 ]','FontSize',20);

%% Plot Vesta

figure; hold on; axis equal

% tri_fig=trisurf(FV,xi,yi,zi,ccol); shading interp
% set(tri_fig,'edgecolor','k','facecolor',ccol(1,:)); 
% trimesh_fig=trimesh(FV,xi,yi,zi,ri,'edgecolor','k','facecolor','none');

p = patch('Faces',FV,'Vertices',[x1 y1 z1],'FaceColor','b');
% p = patch('Faces',FV,'Vertices',[x2 y2 z2],'FaceColor','b');
% p = patch('Faces',FV,'Vertices',[x3 y3 z3],'FaceColor','b');

clear cdata

cdata = rho1_s';
% cdata = (rho1_inc+rho2_inc)';
set(p,'FaceColor','flat',...
'FaceVertexCData',cdata,...
'CDataMapping','scaled')

 caxis([-10000 10000]);

axis equal

cbar=colorbar
ylabel(cbar,'Density [kg/m^3] ','FontSize',20);

% rho(abs((rho-rho_mean))>500)=rho_mean










