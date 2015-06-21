ccc


%  ShapeFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/SHAPE-120709.TXT';
% %ShapeFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/SHAPE-120604.TXT';
% % 
% lmcosi_hamo2=MakeSHShapeModel(ShapeFileName);
% % 
% % 
% % plotplm(lmcosi,[],[],2,.1,[],[],[])
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
N_trunc=16;
[C_obs,S_obs,C_obs_std,S_obs_std]=TruncateGravityModel(C_obs,S_obs,C_obs_std,S_obs_std,N_trunc);
[C_gt,S_gt,C_obs_std,S_obs_std]=TruncateGravityModel(C_gt,S_gt,C_obs_std,S_obs_std,N_trunc);

lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);
lmcosi_gt=TruncateGravityModel(lmcosi_gt,N_trunc,0);

V=7.4969498900204E+16;
ro_mean=3.4546423486191E+03;
G=6.67259e-11;

%% core size
% H
% r_core=113500;
% a_core=117900;
% b_core=a_core;
% c_core=(r_core^3)/(a_core^2);
% f_core=(a_core-c_core)/a_core;


% CR
% r_core=106400;
% a_core=112300;
% b_core=a_core;
% c_core=(r_core^3)/(a_core^2);
% f_core=(a_core-c_core)/a_core;

% CH
% r_core=132300;
% a_core=136600;
% b_core=a_core;
% c_core=(r_core^3)/(a_core^2);
% f_core=(a_core-c_core)/a_core;

% 3/4H + 1/4CM
r_core=114200;
a_core=118807;
b_core=a_core;
c_core=(r_core^3)/(a_core^2);
f_core=(a_core-c_core)/a_core;

% JPL
% r_core=110000;
% a_core=114000;
% b_core=a_core;
% c_core=(r_core^3)/(a_core^2);
% f_core=(a_core-c_core)/a_core;



a_mantle=232217.5;
b_mantle=a_mantle;
c_mantle=199000;

MinDegree=1;
MaxDegree=17;

%% Densities

% H
% ro_crust=2903;
% ro_mantle=3270;
% ro_core=6540;

% CR
% ro_crust=2923;
% ro_mantle=3345;
% ro_core=6175;
% 
% CH
% ro_crust=2642;
% ro_mantle=2985;
% ro_core=7802;

% 3/4H + 1/4CM
ro_crust=2910;
ro_mantle=3294;
ro_core=6259;

% JPL
% ro_crust=2900;
% ro_mantle=3200;
% ro_core=7400;




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

C_core(1,2)=-C_gt(1,2)*mu_crust_diff/mu_core_diff;
S_core(1,2)=-S_gt(1,2)*mu_crust_diff/mu_core_diff;
C_core(1,1)=-C_gt(1,1)*mu_crust_diff/mu_core_diff;

x_core_cm=C_core(1,2)*Rref;
y_core_cm=S_core(1,2)*Rref;
z_core_cm=C_core(1,1)*Rref;

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
Vesta.crust_mass_fraction=1-Vesta.core_mass_fraction-Vesta.mantle_mass_fraction;

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
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,g_ell_v_obs,BouguerFileName);

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

% MapRadialGrid(dg_bouguer_ell_an);

%% Plot mantle and outer shape

fig1=figure;
hold on


Resolution=1;
MaxDegreeTopo=100;

x_mantle=a_mantle_new*cos(lambdai).*cos(fii);
y_mantle=a_mantle_new*sin(lambdai).*cos(fii);
z_mantle=c_mantle_new*sin(fii);

surf(x_mantle,y_mantle,z_mantle);
 
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

plotplm(lmcosi_shape,[],[],2,Resolution,[],[],[]);

axis equal

alpha(0.8);










