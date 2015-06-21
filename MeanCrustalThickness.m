function [mean_crustal_thickness,max_crust,min_crust,max_bouguer,min_bouguer,f_mantle_new]=MeanCrustalThickness(ro_crust,ro_mantle,ro_core,C_obs,S_obs,mu_obs,Rref,lmcosi_obs_full,C_gt,S_gt,lmcosi_gt)
%  ShapeFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/SHAPE-120709.TXT';
% %ShapeFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/SHAPE-120604.TXT';
% % 
% lmcosi_hamo2=MakeSHShapeModel(ShapeFileName);
% % 
% % 
% % plotplm(lmcosi,[],[],2,.1,[],[],[])
% 
% 
% 
% 
%   WriteInputTopoSHFile(lmcosi_hamo2);
% offset=[-0.33169D+00   -0.14037D+01    0.10498D+00];

mu=mu_obs;

%% Truncate spherical harmonic models
N_trunc=16;
[C_obs,S_obs,~,~]=TruncateGravityModel(C_obs,S_obs,C_obs*0,C_obs*0,N_trunc);
[C_gt,S_gt,~,~]=TruncateGravityModel(C_gt,S_gt,C_obs*0,C_obs*0,N_trunc);

lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);
lmcosi_gt=TruncateGravityModel(lmcosi_gt,N_trunc,0);

V=7.4969498900204E+16;
ro_mean=3.4546423486191E+03;
G=6.67259e-11;

a_core=114300;
b_core=a_core;
c_core=101800;

a_mantle=232217.5;
b_mantle=a_mantle;
c_mantle=198000;

MinDegree=1;
MaxDegree=17;

% ro_crust=2900;
% ro_mantle=3200;
% ro_core=7800;

ro_mantle_diff=ro_mantle-ro_crust;
ro_core_diff=ro_core-ro_mantle; 

mu_mantle_diff=muEllipsoid(a_mantle,b_mantle,c_mantle,ro_mantle_diff);
mu_core_diff=muEllipsoid(a_core,b_core,c_core,ro_core_diff);
mu_crust_diff=ro_crust*V*G;

mu_calc=mu_core_diff+mu_mantle_diff+mu_crust_diff;

mu_diff=mu_obs-mu_calc;

mu_mantle_ratio=(mu_mantle_diff+mu_diff)/mu_mantle_diff;

a_mantle=a_mantle*(mu_mantle_ratio^(1/3));
b_mantle=b_mantle*(mu_mantle_ratio^(1/3));
c_mantle=c_mantle*(mu_mantle_ratio^(1/3));

mu_mantle_diff=muEllipsoid(a_mantle,b_mantle,c_mantle,ro_mantle_diff);
mu_core_diff=muEllipsoid(a_core,b_core,c_core,ro_core_diff);
mu_crust_diff=ro_crust*V*G;

mu_calc=mu_core_diff+mu_mantle_diff+mu_crust_diff;

[C_core]=SHRotationalEllipsoid(a_core,c_core,MaxDegree,N_trunc,Rref);
S_core=zeros(N_trunc,N_trunc+1);
S_mantle=zeros(N_trunc,N_trunc+1);

x_cf=C_gt(1,2)*Rref/NormCoef(1,1);
y_cf=S_gt(1,2)*Rref/NormCoef(1,1);
z_cf=C_gt(1,1)*Rref/NormCoef(1,0);


% C_core(1,2)=-C_gt(1,2)*mu_crust_diff/mu_core_diff;
% S_core(1,2)=-S_gt(1,2)*mu_crust_diff/mu_core_diff;
% C_core(1,1)=-C_gt(1,1)*mu_crust_diff/mu_core_diff;
% 
% x_core_cm=C_core(1,2)*Rref;
% y_core_cm=S_core(1,2)*Rref;
% z_core_cm=C_core(1,1)*Rref;


%% Setting core offset
x_core_cm=-830
y_core_cm=-200
z_core_cm=-5660
%


% CgS up -> BAS down -> Crust thickness S up, Crust thickness N downa
% 
C_core(1,2)=x_core_cm/Rref;
S_core(1,2)=y_core_cm/Rref;
C_core(1,1)=z_core_cm/Rref;






disp('Center of figure coordinates')
x_cf
y_cf
z_cf

fi_cf=atan(z_cf/sqrt(x_cf^2+y_cf^2));
lambda_cf=atan2(y_cf,x_cf);

disp('Core center corrdinates')
 x_core_cm
 y_core_cm
 z_core_cm
 
disp('Core offset')
 norm([x_core_cm y_core_cm z_core_cm])

C_20_mantle=(C_obs(2,1)*mu_obs-C_core(2,1)*mu_core_diff-C_gt(2,1)*mu_crust_diff)/mu_mantle_diff;

x=C_20_mantle*5*Rref*Rref/NormCoef(2,0);
y=3*mu_mantle_diff/(4*pi*ro_mantle_diff*G);
[a_mantle_new,c_mantle_new]=FindAxes(x,y);

[C_mantle]=SHRotationalEllipsoid(a_mantle_new,c_mantle_new,MaxDegree,N_trunc,Rref);

% C_mantle(2,1)=C_20_mantle;

f_mantle_new=(a_mantle_new-c_mantle_new)/a_mantle_new;

%% Computing gravity calc
C_calc=(C_core*mu_core_diff+C_mantle*mu_mantle_diff+C_gt*mu_crust_diff)/mu_obs;
S_calc=(S_core*mu_core_diff+S_mantle*mu_mantle_diff+S_gt*mu_crust_diff)/mu_obs;
 
C_calc(2,1)-C_obs(2,1)

Vesta.ro_core=ro_core;
Vesta.ro_mantle=ro_mantle;
Vesta.ro_crust=ro_crust;
Vesta.a_core=a_core;
Vesta.c_core=c_core;
Vesta.a_mantle=a_mantle_new;
Vesta.c_mantle=c_mantle_new;

Vesta.core_mass_fraction=muEllipsoid(a_core,a_core,c_core,ro_core)/mu_obs;
Vesta.mantle_mass_fraction=muEllipsoid(a_mantle_new,a_mantle_new,c_mantle_new,ro_mantle)/mu_obs;
Vesta.crust_fraction=1-Vesta.core_mass_fraction-Vesta.mantle_mass_fraction;

Vesta

ag=2.918299428885863e+05;
bg=2.650067859489697e+05;

FiStep=1;
LambdaStep=1;
[lambdai,fii]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );
x=ag*cos(lambdai).*cos(fii);
y=ag*sin(lambdai).*cos(fii);
z=bg*sin(fii);
r2=x.*x+y.*y+z.*z;

C_bouguer=C_obs-C_calc;
S_bouguer=S_obs-S_calc;

C_bouguer(1:MinDegree-1,:)=0;
S_bouguer(1:MinDegree-1,:)=0;

% dg_bouguer=(RadialGForceFromGravityModel(x,y,z,C_bouguer,S_bouguer,mu_obs,Rref)-mu_obs./r2)*1e5;
% %% Radial Bouguer Anomaly
% BouguerFileName='bouguertest.txt';
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_bouguer,BouguerFileName);
% 
 lmcosi_calc=CS2lmcosi(C_calc,S_calc);
 lmcosi_bouguer=CS2lmcosi(C_bouguer,S_bouguer);

%% Vector Bouguer anomaly
% BouguerFileName='bouguertest_vector.txt';
% 
[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi_obs,x,y,z);
[gx_calc,gy_calc,gz_calc]=GravityAcceleration(mu,Rref,lmcosi_calc,x,y,z);
% 
% g_obs=1e5*sqrt(gx_obs.*gx_obs+gy_obs.*gy_obs+gz_obs.*gz_obs);
% 
% dgx_bouguer=1e5*(gx_obs-gx_calc);
% dgy_bouguer=1e5*(gy_obs-gy_calc);
% dgz_bouguer=1e5*(gz_obs-gz_calc);
% 
% dg_bouguer_vector=sqrt(dgx_bouguer.*dgx_bouguer+dgy_bouguer.*dgy_bouguer+dgz_bouguer.*dgz_bouguer);
% 
% [g_ell_v_obs]=NormalGravityComponent(dgx_bouguer,dgy_bouguer,dgz_bouguer,x,y,z);

% WriteXYZ(lambdai*180/pi,fii*180/pi,g_ell_v_obs,BouguerFileName);

%% Ellipsoidal Bouguer Anomaly
% BouguerFileName='bouguertest_ell.txt';
% 
% [g_ell_obs]=1e5*NormalGravityComponent(gx_obs,gy_obs,gz_obs,x,y,z);
% [g_ell_calc]=1e5*NormalGravityComponent(gx_calc,gy_calc,gz_calc,x,y,z);
% 
% dg_ell_bouguer=(g_ell_obs-g_ell_calc);
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_ell_bouguer,BouguerFileName);

% SphericalHarmonicCorrelation(lmcosi_calc,lmcosi_obs);
% SphericalHarmonicCorrelation(lmcosi_bouguer,lmcosi_obs);

% %%  Bouguer Anomaly Difference
% BouguerFileName='bouguertest_ell_rad_diff.txt';
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_ell_bouguer-dg_bouguer,BouguerFileName);
% 
% BouguerFileName='bouguertest_vec_rad_diff.txt';
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_bouguer_vector-abs(dg_bouguer),BouguerFileName);
% 
% BouguerFileName='bouguertest_vec_ell_diff.txt';
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_bouguer_vector-abs(dg_ell_bouguer),BouguerFileName);

%% Bouguer anomaly analytical ellipsoidal
BouguerFileName='bouguertest_ell_an.txt';

[g_up_obs,g_east_obs,g_north_obs]=GravityComponents(gx_obs,gy_obs,gz_obs,x,y,z,ag,bg);
[g_up_calc,g_east_calc,g_north_calc]=GravityComponents(gx_calc,gy_calc,gz_calc,x,y,z,ag,bg);

g_h_obs=sqrt(g_north_obs.*g_north_obs+g_east_obs.*g_east_obs);

g_up_obs=1e5*g_up_obs;
g_east_obs=1e5*g_east_obs;
g_north_obs=1e5*g_north_obs;

g_up_calc=1e5*g_up_calc;
g_east_calc=1e5*g_east_calc;
g_north_calc=1e5*g_north_calc;

dg_bouguer_ell_an=g_up_obs-g_up_calc;
% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_bouguer_ell_an,BouguerFileName);

%% Map Bouguer anomaly

%  MapRadialGrid(dg_bouguer_ell_an);

%% MinMax Bouguer anomaly
% max_bouguer=max(dg_bouguer_ell_an(:));
% min_bouguer=min(dg_bouguer_ell_an(:));

%% Solving for the mantle interface

eps_ri_mantle=50;
eps_mantle_fa_correction=50;
lambda_lagrange=3.5e-15;
DegreeCrit=5;


% 15.5 ->3
% 3.5 ->5
% 0.38 -> 9
% MaxTopoPower=10;
% MaxIter=10;

[xi_ell_mantle,yi_ell_mantle,zi_ell_mantle]=MakeRotationalEllipsoid(a_mantle_new,c_mantle_new,FiStep,LambdaStep);
[xi_ell_core,yi_ell_core,zi_ell_core]=MakeRotationalEllipsoid(a_core,c_core,FiStep,LambdaStep);

xi_ell_core=xi_ell_core+x_core_cm;
yi_ell_core=yi_ell_core+y_core_cm;
zi_ell_core=zi_ell_core+z_core_cm;

[~,~,ri_mantle]=cart2sph(xi_ell_mantle,yi_ell_mantle,zi_ell_mantle);
[~,~,ri_core]=cart2sph(xi_ell_core,yi_ell_core,zi_ell_core);

[lmcosi_mantle_shape,~]=xyz2plm(ri_mantle,N_trunc,'im',[],[],[]);
[lmcosi_core_shape,~]=xyz2plm(ri_core,N_trunc,'im',[],[],[]);

% mantle_shape_fig=figure('Position',[1 1 1000 1000]);
% plotplm(lmcosi_mantle_shape,[],[],2,1,[],[],[]);

D_mantle=lmcosi_mantle_shape(1,3);

lmcosi_mantle_shape_correction_1=CreateEmptylmcosi(N_trunc);
lmcosi_mantle_shape_fa_correction=CreateEmptylmcosi(N_trunc);
lmcosi_mantle_shape_power=CreateEmptylmcosi(N_trunc);
lmcosi_mantle_shape_new=CreateEmptylmcosi(N_trunc);
lmcosi_bouguer=AddZeroHarm(lmcosi_bouguer,0);

[sdl_bouguer,l_bouguer,~,~,~,~]=plm2spec(lmcosi_bouguer);

% [bouguer_anomaly_fig,bouguer_anomaly_axes]=BlackFigureLog(' ','Power per degree','Degree');
% [sdl_grav,l_grav,~,~,~,~]=plm2spec(lmcosi_obs_full);
% semilogy(l_grav,sdl_grav,'yo-','LineWidth',2);

% [sdl_grav_error,l_grav_error,~,~,~,~]=plm2spec([lmcosi_obs_full(:,1:2) lmcosi_obs_full(:,5:6)]);
% semilogy(l_grav_error,sdl_grav_error,'co-','LineWidth',2);


% semilogy(l_bouguer,sdl_bouguer,'wo-','LineWidth',2);
% hold on;

Degree=lmcosi_mantle_shape(:,1);
M=mu_obs/G;

Coef1=M*(2.*Degree+1).*((Rref./D_mantle).^Degree)./(4*pi*ro_mantle_diff.*(D_mantle^2));


lambda_lagrange=1./(M*(2.*DegreeCrit+1).*((Rref./D_mantle).^DegreeCrit)./(4*pi*ro_mantle_diff.*(D_mantle^2))).^2;


w=(1+lambda_lagrange.*(Coef1.^2)).^(-1);

% [w_fig,w_axes]=BlackFigure('Downward continuation filter','w','Degree');
% w_fig=figure;
% set(gcf, 'Units','centimeters', 'Position',[0 0 13 9]);
% set(gcf, 'PaperPositionMode','auto')
% set(gca,'XTick',0:2:20);
% plot(Degree,w,'g','LineWidth',2);
% xlabel('Degree','FontSize',12);
% ylabel('w','fontsize',12);
% grid on
% 
% print(w_fig, '-dpsc', 'w_filter.eps');

% PrintBlack('w.eps')


lmcosi_mantle_shape_correction_1(:,3:4)=[w,w].*(lmcosi_bouguer(:,3:4)).*[Coef1,Coef1];

lmcosi_mantle_shape_correction_full=lmcosi_mantle_shape_correction_1;

dri_mantle=eps_ri_mantle+1;

iter=1;


while (dri_mantle>eps_ri_mantle)
    

[ri_mantle_shape_correction,~]=plm2xyz(lmcosi_mantle_shape_correction_full,1);

% ri_mantle_topo=ri_mantle-D_mantle;

n=2;

dri_mantle_fa_correction=eps_mantle_fa_correction+1;

while (dri_mantle_fa_correction>eps_mantle_fa_correction)
    
P=ones(numel(Degree),1);       
    
    for j=1:n        
        P=P.*(Degree+4-j);       
    end
        
Coef2=P./((D_mantle^n)*factorial(n).*(Degree+3));
    
ri_mantle_shape_correction_n=ri_mantle_shape_correction.^n;    

[lmcosi_mantle_shape_correction_power,dw]=xyz2plm(ri_mantle_shape_correction_n,N_trunc,'im',[],[],[]);  

lmcosi_mantle_shape_fa_correction(:,3:4)=[w, w].*D_mantle.*lmcosi_mantle_shape_correction_power(:,3:4).*[Coef2, Coef2];
    
lmcosi_mantle_shape_correction_full(:,3:4)=lmcosi_mantle_shape_correction_1(:,3:4)-lmcosi_mantle_shape_fa_correction(:,3:4);

% lmcosi_shape_mantle_new(:,3:4)=lmcosi_shape_mantle(:,3:4)+lmcosi_shape_correction(:,3:4);   
%
% plotplm(lmcosi_shape_mantle_new,[],[],2,1,[],[],[]);
[ri_mantle_fa_correction,~]=plm2xyz(lmcosi_mantle_shape_fa_correction,1);

dri_mantle_fa_correction=max(max(abs(ri_mantle_fa_correction)));

n=n+1;

end

lmcosi_mantle_shape_new(:,3:4)=lmcosi_mantle_shape(:,3:4)+lmcosi_mantle_shape_correction_full(:,3:4);

% plotplm(lmcosi_mantle_shape_new,[],[],2,1,[],[],[]);
[ri_mantle_new,~]=plm2xyz(lmcosi_mantle_shape_new,1);

dri_mantle=max(max(abs(ri_mantle_new-ri_mantle)));

ri_mantle=ri_mantle_new;

iter=iter+1;

end

%% Check new bouguer anomaly
MaxTopoPower=15;
MaxDegreeTopo=16;
MaxDegreeGrav=16;

[ri_mantle_new,~,~,~]=plm2xyz(lmcosi_mantle_shape_new,1);

lmcosi_mantle_grav_new=TopoSH2GravitySH(ri_mantle_new,mu_mantle_diff,ro_mantle_diff,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_mantle_grav_new(1,:)=[];

[C_mantle_corr,S_mantle_corr]=lmcosi2CS(lmcosi_mantle_grav_new);

C_calc_new=(C_core*mu_core_diff+C_mantle_corr*mu_mantle_diff+C_gt*mu_crust_diff)/mu_obs;
S_calc_new=(S_core*mu_core_diff+S_mantle_corr*mu_mantle_diff+S_gt*mu_crust_diff)/mu_obs;

lmcosi_calc_new=CS2lmcosi(C_calc_new,S_calc_new);

[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi_obs,x,y,z);
[gx_calc_new,gy_calc_new,gz_calc_new]=GravityAcceleration(mu,Rref,lmcosi_calc_new,x,y,z);

BouguerFileName='bouguertest_ell_corr.txt';

[g_up_calc_new,g_east_calc_new,g_north_calc_new]=GravityComponents(gx_calc_new,gy_calc_new,gz_calc_new,x,y,z,ag,bg);
[g_up_obs,g_east_obs,g_north_obs]=GravityComponents(gx_obs,gy_obs,gz_obs,x,y,z,ag,bg);

g_up_calc_new=1e5*g_up_calc_new;
g_up_obs=1e5*g_up_obs;

g_east_calc_new=1e5*g_east_calc_new;
g_north_calc_new=1e5*g_north_calc_new;

dg_bouguer_ell_an_new=g_up_obs-g_up_calc_new;

%% MinMax minimized Bouguer anomaly

max_bouguer=max(dg_bouguer_ell_an_new(:));
min_bouguer=min(dg_bouguer_ell_an_new(:));


% WriteXYZ(lambdai*180/pi,fii*180/pi,dg_bouguer_ell_an_new,BouguerFileName);

% MapRadialGrid(dg_bouguer_ell_an_new);
% 
% close all
 
C_bouguer=C_obs-C_calc_new;
S_bouguer=S_obs-S_calc_new;

C_bouguer(1:MinDegree-1,:)=0;
S_bouguer(1:MinDegree-1,:)=0;

lmcosi_calc=lmcosi_calc_new;
lmcosi_bouguer=CS2lmcosi(C_bouguer,S_bouguer);
lmcosi_bouguer=AddZeroHarm(lmcosi_bouguer,0);


%% Map mantle and shape

% reference ellipsoid
a_ref_ell=280917.27;
c_ref_ell=226248.16;

Resolution=1;
MaxDegreeTopo=100;

% int_str_fig=figure('color',[0 0 0]);
% plotplm(lmcosi_mantle_shape_new,[],[],2,1,[],[],[]);
% alpha(0.5)
% hold on
%  
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
%  lmcosi_shape(4,3)=0;
%  lmcosi_shape(1,3)=0;
% plotplm(lmcosi_shape,[],[],2,Resolution,[],[],[]);
% alpha(0.4)
% hold on

[ri_mantle,~,~,~]=plm2xyz(lmcosi_mantle_shape,Resolution);
[ri_mantle_new,lambda_mantle_new,fi_mantle_new,~]=plm2xyz(lmcosi_mantle_shape_new,Resolution, [0 90 360 -90]);
[ri_shape,lambda_shape,fi_shape,~]=plm2xyz(lmcosi_shape,Resolution, [0 90 360 -90]);


[lambda_mantle_new,fi_mantle_new]=meshgrid(lambda_mantle_new,fi_mantle_new);
[lambda_shape,fi_shape]=meshgrid(lambda_shape,fi_shape);

[x_mantle_new,y_mantle_new,z_mantle_new]=sph2cart(lambda_mantle_new/180*pi,fi_mantle_new/180*pi,ri_mantle_new);
[x_shape,y_shape,z_shape]=sph2cart(lambda_shape/180*pi,fi_shape/180*pi,ri_shape);




[~,~,H_mantle_new]=XYZ2BLH(x_mantle_new,y_mantle_new,z_mantle_new,a_ref_ell,Eccentricity(a_ref_ell,c_ref_ell));
[~,~,H_shape]=XYZ2BLH(x_shape,y_shape,z_shape,a_ref_ell,Eccentricity(a_ref_ell,c_ref_ell));



% MapRadialGrid(flipud(ri_shape-ri_mantle_new)/1000);

% title('Crustal thickness','fontsize',25);
% hold on;
% WriteXYZ(lambdai*180/pi,fii*180/pi,flipud(ri_shape-ri_mantle_new)/1000,'Crustal_Thickness.xyz');

% WriteXYZ(lambdai*180/pi,fii*180/pi,flipud(H_shape-H_mantle_new)/1000,'Crustal_Thickness.xyz');


[fii_rand,lambdai_rand]=GenerateRandomSphericalCoord(2000);

H_shape_rand=griddata(lambda_shape-180,fi_shape,H_shape,lambdai_rand*180/pi,fii_rand*180/pi);
H_mantle_new_rand=griddata(lambda_mantle_new-180,fi_mantle_new,H_mantle_new,lambdai_rand*180/pi,fii_rand*180/pi);

crustal_thickness=(H_shape_rand-H_mantle_new_rand)/1000;

mean_crustal_thickness=mean(crustal_thickness(:));

max_crust=max(crustal_thickness(:));
min_crust=min(crustal_thickness(:));

