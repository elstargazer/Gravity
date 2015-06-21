function [angle_gc,dg_bouguer_ell_an_gc,mean_crustal_thickness]=BouguerAnomalyDH(ro_crust,ro_mantle) 
% close all;
% clc
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

[C_gt,S_gt]=ReadBalminoSH('file_harmo_pot');
lmcosi_gt=CS2lmcosi(C_gt,S_gt);

Rref=265000;
% C_gt(1,2)=offset(1)/Rref;
% S_gt(1,2)=offset(2)/Rref;
% C_gt(1,1)=offset(3)/Rref;

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

[C_obs,S_obs,C_obs_std,S_obs_std,mu_obs,Rref]=LoadGravityModel(GravityFileName);

%% Truncate spherical harmonic models
N_trunc=15;
[C_obs,S_obs,C_obs_std,S_obs_std]=TruncateGravityModel(C_obs,S_obs,C_obs_std,S_obs_std,N_trunc);
[C_gt,S_gt,C_obs_std,S_obs_std]=TruncateGravityModel(C_gt,S_gt,C_obs_std,S_obs_std,N_trunc);

lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);
lmcosi_gt=TruncateGravityModel(lmcosi_gt,N_trunc,0);

V=7.4969498900204E+16;
ro_mean=3.4546423486191E+03;
G=6.67259e-11;

a_core=118807;
b_core=a_core;
c_core=1.055150119504450e+05;

a_mantle=232217.5;
b_mantle=a_mantle;
c_mantle=198000;

MinDegree=1;
MaxDegree=15;

Factor=1;

% ro_crust=2910*Factor;
% ro_mantle=3294*Factor;
ro_core=6259;

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

%% Computing core offset
C_core(1,2)=-C_gt(1,2)*mu_crust_diff/mu_core_diff;
S_core(1,2)=-S_gt(1,2)*mu_crust_diff/mu_core_diff;
C_core(1,1)=-C_gt(1,1)*mu_crust_diff/mu_core_diff;

x_core_cm=C_core(1,2)*Rref;
y_core_cm=S_core(1,2)*Rref;
z_core_cm=C_core(1,1)*Rref;

%% Setting core offset
% x_core_cm=-830;
% y_core_cm=-200;
% z_core_cm=-5660;
% 
% 
% C_core(1,2)=x_core_cm/Rref;
% S_core(1,2)=y_core_cm/Rref;
% C_core(1,1)=z_core_cm/Rref;


%% Printing core offset

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
Vesta.mantle_crust_fraction=1-Vesta.core_mass_fraction-Vesta.mantle_mass_fraction;

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
BouguerFileName='bouguertest_vector.txt';

[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi_obs,x,y,z);
[gx_calc,gy_calc,gz_calc]=GravityAcceleration(mu,Rref,lmcosi_calc,x,y,z);

g_obs=1e5*sqrt(gx_obs.*gx_obs+gy_obs.*gy_obs+gz_obs.*gz_obs);

dgx_bouguer=1e5*(gx_obs-gx_calc);
dgy_bouguer=1e5*(gy_obs-gy_calc);
dgz_bouguer=1e5*(gz_obs-gz_calc);

dg_bouguer_vector=sqrt(dgx_bouguer.*dgx_bouguer+dgy_bouguer.*dgy_bouguer+dgz_bouguer.*dgz_bouguer);

[g_ell_v_obs]=NormalGravityComponent(dgx_bouguer,dgy_bouguer,dgz_bouguer,x,y,z);



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

%% Map Bouguer anomaly

% MapRadialGrid(dg_bouguer_ell_an);

MinBouguerAnom=min(dg_bouguer_ell_an(:));
MaxBouguerAnom=max(dg_bouguer_ell_an(:));

% caxis([MinBouguerAnom MaxBouguerAnom]);

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



Degree=lmcosi_mantle_shape(:,1);
M=mu_obs/G;

Coef1=M*(2.*Degree+1).*((Rref./D_mantle).^Degree)./(4*pi*ro_mantle_diff.*(D_mantle^2));

lambda_lagrange=1./(M*(2.*DegreeCrit+1).*((Rref./D_mantle).^DegreeCrit)./(4*pi*ro_mantle_diff.*(D_mantle^2))).^2;

w=(1+lambda_lagrange.*(Coef1.^2)).^(-1);

% [w_fig,w_axes]=BlackFigure('Downward continuation filter','w','Degree');


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
MaxDegreeTopo=100;
MaxDegreeGrav=15;

[ri_mantle_new,~,~,~]=plm2xyz(lmcosi_mantle_shape_new,1);

lmcosi_mantle_grav_new=TopoSH2GravitySH(ri_mantle_new,mu_mantle_diff,ro_mantle_diff,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_mantle_grav_new(1,:)=[];

[C_mantle_corr,S_mantle_corr]=lmcosi2CS(lmcosi_mantle_grav_new);

C_calc_new=(C_core*mu_core_diff+C_mantle_corr*mu_mantle_diff+C_gt*mu_crust_diff)/mu_obs;
S_calc_new=(S_core*mu_core_diff+S_mantle_corr*mu_mantle_diff+S_gt*mu_crust_diff)/mu_obs;



C_core_per=(C_core*mu_core_diff/mu_obs)./C_calc_new*100
S_core_per=(S_core*mu_core_diff/mu_obs)./S_calc_new*100

C_mantle_per=(C_mantle_corr*mu_mantle_diff/mu_obs)./C_calc_new*100
S_mantle_per=(S_mantle_corr*mu_mantle_diff/mu_obs)./S_calc_new*100


C_shape_per=(C_gt*mu_crust_diff/mu_obs)./C_calc_new*100
S_shape_per=(S_gt*mu_crust_diff/mu_obs)./S_calc_new*100

lmcosi_calc_new=CS2lmcosi(C_calc_new,S_calc_new);
lmcosi_core_grav_new=CS2lmcosi(C_core,S_core);


lmcosi_calc_new2=TruncateGravityModel(lmcosi_calc_new,15,0);

[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi_obs,x,y,z);
[gx_calc_new,gy_calc_new,gz_calc_new]=GravityAcceleration(mu,Rref,lmcosi_calc_new,x,y,z);



[g_up_calc_new,g_east_calc_new,g_north_calc_new]=GravityComponents(gx_calc_new,gy_calc_new,gz_calc_new,x,y,z,ag,bg);

g_up_calc_new=1e5*g_up_calc_new;
g_east_calc_new=1e5*g_east_calc_new;
g_north_calc_new=1e5*g_north_calc_new;

dg_bouguer_ell_an_new=g_up_obs-g_up_calc_new;

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

% MapRadialGrid((ri_shape));


% WriteXYZ(lambdai*180/pi,fii*180/pi,flipud(ri_shape-ri_mantle_new)/1000,'Crustal_Thickness.xyz');



[fii_rand,lambdai_rand]=GenerateRandomSphericalCoord(2000);

H_shape_rand=griddata(lambda_shape-180,fi_shape,H_shape,lambdai_rand*180/pi,fii_rand*180/pi);
H_mantle_new_rand=griddata(lambda_mantle_new-180,fi_mantle_new,H_mantle_new,lambdai_rand*180/pi,fii_rand*180/pi);

mean_crustal_thickness=mean(H_shape_rand(:)-H_mantle_new_rand(:))/1000




%% Profile

fi_rs=-60.0251;
lambda_rs=-50;

az1=180;
az2=az1+180;
rng=90;
Npts=100;

[fi_gc1,lambda_gc1] = track1(fi_rs,lambda_rs,az1,rng);
[fi_gc2,lambda_gc2] = track1(fi_rs,lambda_rs,az2,rng);



% plotm(fi_rs/180*pi,lambda_rs/180*pi,'or','MarkerSize',10);
% 
% plotm(fi_gc1/180*pi,lambda_gc1/180*pi,'k-','LineWidth',3);
% plotm(fi_gc2/180*pi,lambda_gc2/180*pi,'k-','LineWidth',3);

fi_gc=[flipud(fi_gc1); fi_gc2];
lambda_gc=[flipud(lambda_gc1); lambda_gc2];


angle_gc=linspace(-90,90,2*Npts)';

% [fig_profile,fig_profile_ax]=BlackFigure(' ','Ellipsoidal height','Degrees along the profile');




%% Bouguer Profile
Coef=1.00;
x_ell_gc=Coef*ag*cos(lambda_gc/180*pi).*cos(fi_gc/180*pi);
y_ell_gc=Coef*ag*sin(lambda_gc/180*pi).*cos(fi_gc/180*pi);
z_ell_gc=Coef*bg*sin(fi_gc/180*pi);

[gx_obs_gc,gy_obs_gc,gz_obs_gc]=GravityAcceleration(mu,Rref,lmcosi_obs,x_ell_gc,y_ell_gc,z_ell_gc);
[gx_calc_gc,gy_calc_gc,gz_calc_gc]=GravityAcceleration(mu,Rref,lmcosi_calc,x_ell_gc,y_ell_gc,z_ell_gc);

[g_up_obs_gc,g_east_obs_gc,g_north_obs_gc]=GravityComponents(gx_obs_gc,gy_obs_gc,gz_obs_gc,x_ell_gc,y_ell_gc,z_ell_gc,Coef*ag,Coef*bg);
[g_up_calc_gc,g_east_calc_gc,g_north_calc_gc]=GravityComponents(gx_calc_gc,gy_calc_gc,gz_calc_gc,x_ell_gc,y_ell_gc,z_ell_gc,Coef*ag,Coef*bg);

dg_bouguer_ell_an_gc=(g_up_obs_gc-g_up_calc_gc);

% [fig_bouguer_profile,fig_bouguer_profile_ax]=BlackFigure(' ','Bouguer Anomaly [mGal]','Degrees along the profile');

% fig_bouguer_profile=figure;
% set(gca,'FontSize',12);
% 
% hold on;
% plot(angle_gc,dg_bouguer_ell_an_gc*1e5,'b-','LineWidth',3);
% % plot(angle_gc_iv,bouguer_anomaly_gc*1e5,'r-','LineWidth',3);
% 
% % plot(angle_gc, ,'r-','LineWidth',3);
% set(gca,'XTick',-90:15:90);
% axis([-90 90 -200 200]);
% 
% % set(gca, 'Position',[0 0 1 1])
% set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
% set(gcf, 'PaperPositionMode','auto')
% 
% grid off;
% 
% box on;
% 
% xlabel('Degrees along the profile','FontSize',12);
% ylabel('Bouguer anomaly [mGal]','FontSize',12);
% 
% print(fig_bouguer_profile, '-dpsc', 'Bouguer_Rheasilvia_Profile.eps');






