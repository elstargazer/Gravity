ccc


%% Moon gravity

GravityFileName='/Users/antonermakov/GRAIL/Gravity_Models/JGGRAIL_900C9A_SHA.TAB.txt';
% GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

MaxDegreeTopo=100;
MaxDegreeGrav=100;
% Rref_gravity=293000;
Rref_gravity=1738000;
rhomean=3344;
Resolution=1;
MaxTopoPower=7;
L=20;
MinConcentration=0.95;
NTess=2;
circle_rad=30;
MaxDegreeExp=7;
% rho_mean=3457.5;
GoodDegrees=2+L:MaxDegreeGrav-L;

%% Load gravity model

[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);


%% Moon topography
MoonFileName='/Users/antonermakov/GRAIL/Topography/SH/LRO_LTM05_2050_SHA.TAB.txt';
lmcosi_Moon_shape=load(MoonFileName);


q_topo=IsotropicRatio(lmcosi_Moon_shape(2:end,:),lmcosi_Moon_shape(2:end,:));
q_grav=IsotropicRatio(lmcosi_grav,lmcosi_grav);
q_cross=IsotropicRatio(lmcosi_Moon_shape(2:end,:),lmcosi_grav);

figure; hold on;
plot(1:numel(q_topo),q_topo,'-b');
plot(1:numel(q_grav),q_grav,'-r');
plot(1:numel(q_cross),q_cross,'-g');

hold on;
set(gca,'FontSize',20);
xlabel('Degree','FontSize',12);
ylabel('Isotropic ratio','FontSize',20);
xlabel('Degree','FontSize',20);
legend({'Topography','Gravity','Cross'});
ylim([0 2]);
box on;

%% Vesta shape

load Shape1800.mat

q_g=IsotropicRatio(lmcosi_g(2:end,:),lmcosi_g(2:end,:));
q_d=IsotropicRatio(lmcosi_d(2:end,:),lmcosi_d(2:end,:));

figure; hold on;
plot(1:numel(q_d),q_d,'-b');
plot(1:numel(q_g),q_g,'-r');

hold on;
set(gca,'FontSize',20);
xlabel('Degree','FontSize',12);
ylabel('Isotropic ratio','FontSize',20);
xlabel('Degree','FontSize',20);
ylim([0 2]);
legend({'DLR','Gaskell'},'FontSize',20)
box on;


%% Vesta gravity

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);
q=IsotropicRatio(lmcosi_obs_full,lmcosi_obs_full);
figure;
plot(1:numel(q),q,'-b')
hold on
set(gca,'FontSize',20);
xlabel('Degree','FontSize',12);
ylabel('Isotropic ratio','FontSize',20);
xlabel('Degree','FontSize',20);
% ylim([0 2])