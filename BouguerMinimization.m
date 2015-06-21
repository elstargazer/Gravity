ccc
Rref=265000;

V=7.4969498900204E+16;
ro_mean=3.4546423486191E+03;
G=6.67259e-11;

a_core=115000;
b_core=a_core;
c_core=105000;

a_mantle=232217.5;
b_mantle=a_mantle;
c_mantle=198000;

MinDegree=1;
MaxDegree=16;

ro_crust=2990;
ro_mantle=3170;
ro_core=7400;

N_trunc=16;

%% Read gravity from topography

lmcosi_gt=ReadBalminoSH('file_harmo_pot');

lmcosi_gt=TruncateGravityModel(lmcosi_gt,N_trunc);


%% Read observed gravity

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20E/JGV20E02.SHA';
[lmcosi_obs,Rref,mu_obs,mu_obs_std]=ReadGRAILGravityModel(GravityFileName);

lmcosi_obs=TruncateGravityModel(lmcosi_obs,N_trunc);

%% Compute mu's

ro_mantle_diff=ro_mantle-ro_crust;
ro_core_diff=ro_core-ro_mantle; 

mu_mantle_diff=muEllipsoid(a_mantle,b_mantle,c_mantle,ro_mantle_diff);
mu_core_diff=muEllipsoid(a_core,b_core,c_core,ro_core_diff);
mu_crust_diff=ro_crust*V*G;

%% Compute mantle gravity


lmcosi_mantle=lmcosiRotationalEllipsoid(a_mantle,c_mantle,N_trunc,N_trunc,Rref);


lmcosi_core=lmcosiRotationalEllipsoid(a_core,c_core,N_trunc,N_trunc,Rref);


%% Compute gravity from model

lmcosi_calc(:,1:2)=lmcosi_obs(:,1:2);
mu_calc=mu_core_diff+mu_mantle_diff+mu_crust_diff;

lmcosi_calc(:,3:4)=(lmcosi_core(:,3:4)*mu_core_diff+lmcosi_mantle(:,3:4)*mu_mantle_diff+lmcosi_gt(:,3:4)*mu_crust_diff)/...
    (mu_calc);


%% Compute Bouguer coefficients

mu_bouguer=mu_obs-mu_calc;





