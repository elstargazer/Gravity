ccc
tic

%% Physical parameters 

omega=3.267104402965269e-04; % rotation rate of Vesta

%% Search Parameters

N_trunc=16;
MaxDerOrder=7;
eps=200;
Resolution=1;

%% Load gravity model
GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

[lmcosi,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

%% Load topography model

load('VestaHASTALAVESTAshape_sh20.mat');
[ri,lambdai,fii,Plm]=plm2xyz(lmcosi_shape,Resolution,[],[],[]);
[lambdai,fii]=meshgrid(lambdai,fii);
lambdai=lambdai/180*pi;
fii=fii/180*pi;


% FiStep=.2;
% LambdaStep=.2;
% filename='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/HASTA_LA_VESTA_SHAPE/SHAPE-120709.TXT';
% [ri,fii,lambdai]=Load3ColumnsXYZShape(filename,FiStep,LambdaStep);
% 
% [xi,yi,zi]=sph2cart(lambdai,fii,ri);


%% Truncate gravity model

lmcosi=TruncateGravityModel(lmcosi,N_trunc);

lmcosi=AddZeroHarm(lmcosi);

r_b=max(max(ri));

%% Compute potential on shape

U_rot=0.5*omega*omega*(xi.^2+yi.^2);
U_grav=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,xi,yi,zi);
U=U_grav+U_rot;

U_mean=mean(mean(U))

dU=U-U_mean;
gamma=U./ri;
dr=dU./gamma;
ri_g=ri+dr;


%% Plot Vesta
fig=figure('Position',[1 1 1000 1000]);
shape_fig=surf(xi,yi,zi);
set(shape_fig,'FaceColor','none','EdgeColor','k')
axis equal
hold on;

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
    U_grav=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,xi_g,yi_g,zi_g);
    U=U_grav+U_rot;
    dU=U-U_mean;
    gamma=U./ri_g;
    dr=dU./gamma;
    ri_g=ri_g+dr;
    
    iter=iter+1;
    disp(['max dr (' num2str(iter) ')']);
    max(max(abs(dr)))
%     unplot;
    
end


%% Plot geoid
geoid_plot=surf(xi_g,yi_g,zi_g,ri_g);
shading interp
lighting phong
axis equal
alpha(0.4);
 
WriteXYZ(lambdai*180/pi,fii*180/pi,(ri-ri_g)/1000,'HeightWRTGeoid.txt');

%% Find ellipsoid
ell_approx=[281000 226000]


[ell_g,ell_g_rms]=FindRotationalEllipsoid(ri_g,fii,lambdai,ell_approx)
[ell,ell_rms]=FindRotationalEllipsoid(ri,fii,lambdai,ell_approx)



%% Geoid in spherical harmonics

% 
% MaxDegree=20;
% 
% [lmcosi_geoid,dw]=xyz2plm(ri_g,MaxDegree,'im',[],[],[]);

     
% AGUaxes
% surfm(fii,lambdai,ri-ri_g);
% shading interp
% lighting phong

toc 





    
    
    
    