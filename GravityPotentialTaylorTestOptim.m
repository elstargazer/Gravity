% function [fii,lambdai,ri_g,rms]=GravityPotentialTaylorTestOptim(Coef)
ccc
Coef=1;
tic

%% Physical parameters 

omega=3.267104402965269e-04; % rotation rate of Vesta

%% Search Parameters

N_trunc=16;
MaxDerOrder=7;
eps=5;
Npoints=10000;
Resolution=1;
MaxDegree=15;


%% Load gravity model
GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

[lmcosi,Rref,mu,~]=ReadGRAILGravityModel(GravityFileName);

%% Load topography model


load VestaHASTALAVESTAshape_sh720.mat
MaxDegreeTopo=80;
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);



[ri,lambdai,fii,~]=plm2xyz(lmcosi_shape,Resolution,[],[],[]);

[lambdai,fii]=meshgrid(lambdai,fii);

lambdai=lambdai/180*pi;
fii=fii/180*pi;

[xi,yi,zi]=sph2cart(lambdai,fii,ri);

%% General random directions

[fi_rand,lambda_rand]=GenerateRandomSphericalCoord(Npoints);

[r_rand,~,~,~]=plm2xyz(lmcosi_shape,fi_rand*180/pi,lambda_rand*180/pi,[],[]);

[x_rand,y_rand,z_rand]=sph2cart(lambda_rand,fi_rand,r_rand);


% plot3(x_rand,y_rand,z_rand,'.k','MarkerSize',2);
% 
% axis equal
%% Truncate gravity model

lmcosi=TruncateGravityModel(lmcosi,N_trunc,0);

lmcosi=AddZeroHarm(lmcosi,1);

r_b=max(max(ri));


%% Compute potential on shape

U_rot=0.5*omega*omega*(x_rand.^2+y_rand.^2);
U_grav=GravityPotential(mu,Rref,lmcosi,x_rand,y_rand,z_rand);
% U_grav=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,x_rand,y_rand,z_rand);
U=U_grav+U_rot;

U_mean=Coef*mean(mean(U));

U_rot=0.5*omega*omega*(xi.^2+yi.^2);

U_grav=GravityPotential(mu,Rref,lmcosi,xi,yi,zi);
% U_grav=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,xi,yi,zi);
U=U_grav+U_rot;

dU=U-U_mean;
gamma=-U./ri;
dr=-dU./gamma;
ri_g=ri+dr;

%% Plot Vesta
% fig=figure('Position',[1 1 1000 1000]);
% shape_fig=surf(xi,yi,zi);
% set(shape_fig,'FaceColor','none','EdgeColor','k')
% axis equal
% hold on;

%% Compute equipotential sufrace

iter=1;

while (max(max(abs(dr)))>eps)
    
    [xi_g,yi_g,zi_g]=sph2cart(lambdai,fii,ri_g);
    
%     geoid_plot=surf(xi_g,yi_g,zi_g,ri_g);
%     shading interp
%     lighting phong
%     axis equal    
%     alpha(0.4);
    
    U_rot=0.5*omega*omega*(xi_g.^2+yi_g.^2);
    
    U_grav=GravityPotential(mu,Rref,lmcosi,xi_g,yi_g,zi_g);
%     U_grav=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,xi_g,yi_g,zi_g);
    U=U_grav+U_rot;
    dU=U-U_mean;
    gamma=-U./ri_g;
    dr=-dU./gamma;
    ri_g=ri_g+dr;
    
    iter=iter+1;
    disp(['max dr (' num2str(iter) ')']);
    max(max(abs(dr)))
%     unplot;
    
end

%% Plot geoid
% geoid_plot=surf(xi_g,yi_g,zi_g,ri_g);
% shading interp
% lighting phong
% axis equal
% alpha(0.4);
 
WriteXYZ(lambdai*180/pi,fii*180/pi,(ri-ri_g)/1000,'HeightWRTGeoid.txt');

%% Find ellipsoid

ell_approx=[279000 231000];
ell_approx_tri=[276000 276000 231000];
 
% ell_g=[ 276843.07     231464.81]
[ell_g,ell_g_rms]=FindRotationalEllipsoid(ri_g,fii,lambdai,ell_approx)

[ell_g_tri,rms_tri]=FindTriaxialEllipsoidSpice(ri_g,fii,lambdai,ell_approx_tri)
% ell_g_tri = [2.773881791265943   2.762843616566498   2.314811245738660]

[ell,ell_rms]=FindRotationalEllipsoid(ri,fii,lambdai,ell_approx)

ag=ell_g(1);
cg=ell_g(2);

fg=(ag-cg)/ag;

at=ell(1);
ct=ell(2);

ft=(at-ct)/at;

[ft fg]

%% Geoid in spherical harmonics
 
[lmcosi_geoid,dw]=xyz2plm(ri_g,MaxDegree,'im',[],[],[]);
[r_g_rand,~,~,~]=plm2xyz(lmcosi_geoid,fi_rand*180/pi,lambda_rand*180/pi,[],[]);

plotplm(lmcosi_geoid,[],[],2,1,[],[],[]);

[ri,lambdai,fii,~]=plm2xyz(lmcosi_shape,0.1,[],[],[]);
[lambdai,fii]=meshgrid(lambdai,fii);

[ri_geoid,lambdai_geoid,fii_geoid]=plm2xyz(lmcosi_geoid,0.1);
[lambdai_geoid,fii_geoid]=meshgrid(lambdai_geoid,fii_geoid);

[xi_geoid,yi_geoid,zi_geoid]=sph2cart(lambdai_geoid,fii_geoid,ri_geoid);
[xi,yi,zi]=sph2cart(lambdai,fii,ri);

[~,~,H_geoid]=XYZ2BLH(xi_geoid,yi_geoid,zi_geoid,ag,Eccentricity(ag,cg));
[~,~,H]=XYZ2BLH(xi,yi,zi,ag,Eccentricity(ag,cg));

Height_Geoidal=H-H_geoid;

WriteXYZ(lambdai,fii,Height_Geoidal/1000,'VestaHASTALAVESTAshape_geoidal_height.xyz')

 %% rms

rms=sum(sum((r_rand-r_g_rand).^2));
 
 
% AGUaxes
% surfm(fii,lambdai,ri-ri_g);
% shading interp
% lighting phong

toc





    
    
    
    