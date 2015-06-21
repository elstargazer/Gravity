ccc


%% Admittance


GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

lmcosi_grav=AddZeroHarm(lmcosi_grav,1);

Resolution=1;

MaxTopoPower=10;
MaxDegreeTopo=100;
MaxDegreeGrav=20;

% int_str_fig=figure('color',[0 0 0]);
% plotplm(lmcosi_mantle_shape_new,[],[],2,1,[],[],[]);
% alpha(0.5)
% hold on
%  
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

MeanRadius=lmcosi_shape(1,3);

% lmcosi_shape(:,3)=lmcosi_shape(:,3)/MeanRadius;
% lmcosi_shape(:,4)=lmcosi_shape(:,4)/MeanRadius;

%% Computing admittance


% Adm=SphericalHarmonicAdmittance(lmcosi_grav,lmcosi_shape,'r');
% 
% hold on;
% 
% n=1:MaxDegreeTopo;
% 
% CompDepth=10000;
% 
% zeta=CompDepth/MeanRadius;
% 
% rho_ratio=3000/3457;
% 
% AdmTheor=3./(2*n+1).*((1-zeta).^n).*rho_ratio;
% 
% 
%ri=plm2xyz(lmcosi_shape,1)-MeanRadius;
% 
%ri_power=ri.^2;
% 
%lmcosi_shape_power=xyz2plm(ri_power,MaxDegreeTopo,'im',[],[],[]);  
%  
% 
% 
% plot(n,AdmTheor);
% 
% set(gca,'YScale','Log');
% set(gca,'XScale','Log');
% 
% xlim([2 20]);
% ylim([0 0.6]);

%% Admittance of a non-spherical body

Rref=265000;
mu=17.29e9;
rho=3457;


lmcosi_shape(2:end,3)=lmcosi_shape(2:end,3);
lmcosi_shape(2:end,4)=lmcosi_shape(2:end,4);

[ri,~,~,~]=plm2xyz(lmcosi_shape,Resolution);

lmcosi_shape(:,3)=lmcosi_shape(:,3)/1000;
lmcosi_shape(:,4)=lmcosi_shape(:,4)/1000;

lmcosi_shape(:,3)=lmcosi_shape(:,3);
lmcosi_shape(:,4)=lmcosi_shape(:,4);

lmcosi_grav_shape=TopoSH2GravitySH(ri,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

% lmcosi_calc=AddZeroHarm(lmcosi_calc,1);

% 
lmcosi_grav_shape(:,3)=lmcosi_grav_shape(:,3).*(lmcosi_grav_shape(:,1)+1)*mu/(Rref^2)*1e5;
lmcosi_grav_shape(:,4)=lmcosi_grav_shape(:,4).*(lmcosi_grav_shape(:,1)+1)*mu/(Rref^2)*1e5;

lmcosi_grav(:,3)=lmcosi_grav(:,3).*(lmcosi_grav(:,1)+1)*mu/(Rref^2)*1e5;
lmcosi_grav(:,4)=lmcosi_grav(:,4).*(lmcosi_grav(:,1)+1)*mu/(Rref^2)*1e5;
% 
% lmcosi_calc(:,3)=lmcosi_calc(:,3).*(lmcosi_calc(:,1)+1)*mu/(Rref^2)*1e5;
% lmcosi_calc(:,4)=lmcosi_calc(:,4).*(lmcosi_calc(:,1)+1)*mu/(Rref^2)*1e5;

Adm_homo=SphericalHarmonicAdmittance(lmcosi_grav_shape,lmcosi_shape);
Adm_obs=SphericalHarmonicAdmittance(lmcosi_grav,lmcosi_shape);


% Adm_calc=SphericalHarmonicAdmittance(lmcosi_calc,lmcosi_shape);


% set(gca,'YScale','Log');
% set(gca,'XScale','Log');


n=1:MaxDegreeGrav;

hold on;

plot(n,Adm_homo,'o-','Color','r','LineWidth',4,'MarkerSize',2);
plot(n,Adm_obs,'o-','Color','b','LineWidth',4,'MarkerSize',2);



title('Gravity/Topography Admittance','FontSize',25,'FontName','Times');
xlabel('Degree','FontSize',25,'FontName','Times');

grid on;
hold on;

AdmTheor=3./(2*(n)+1).*((MeanRadius/Rref).^n).*(n+1)*mu/(Rref^2)*1e5/1000;

plot(n,AdmTheor,'Color','k','LineWidth',4,'MarkerSize',2);



% plot(1:16,Adm_calc,'o-','Color','g','LineWidth',4,'MarkerSize',2);


% ylim([0 1]);

%%

legend({'homogenous','observed','linear','model'},'fontsize',20);


set(gca,'FontSize',20)

% set(gca,'YScale','Log');

% ylim([0.04,1]);

