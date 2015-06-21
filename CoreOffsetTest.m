


ro_corei=6500:100:8000;

for i=1:numel(ro_corei)
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
N_trunc=16;
[C_obs,S_obs,C_obs_std,S_obs_std]=TruncateGravityModel(C_obs,S_obs,C_obs_std,S_obs_std,N_trunc);
[C_gt,S_gt,C_obs_std,S_obs_std]=TruncateGravityModel(C_gt,S_gt,C_obs_std,S_obs_std,N_trunc);

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

ro_crust=2900;
ro_mantle=3200;
ro_core=ro_corei(i);

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

x_core_cm(i)=C_core(1,2)*Rref;
y_core_cm(i)=S_core(1,2)*Rref;
z_core_cm(i)=C_core(1,1)*Rref;

disp('Center of figure coordinates')
x_cf
y_cf
z_cf

fi_cf=atan(z_cf/sqrt(x_cf^2+y_cf^2));
lambda_cf=atan2(y_cf,x_cf);

end


plot_core_cm = figure('XVisual',''); 
axes1 = axes('Parent',plot_core_cm,'FontSize',12);

hold on

box(axes1,'on');
plot(ro_corei,x_core_cm/1000,'-r','LineWidth',2);
plot(ro_corei,y_core_cm/1000,'-g','LineWidth',2);
plot(ro_corei,z_core_cm/1000,'-b','LineWidth',2);

xlim([min(ro_corei) max(ro_corei)])
ylim([-.1 11.4])

plot(ro_corei,sqrt(x_core_cm.*x_core_cm+y_core_cm.*y_core_cm+z_core_cm.*z_core_cm)/1000,'-k','LineWidth',2);

% set(gca, 'Position',[0 0 1 1])
set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

xlabel('Core density [kg/m^3]','FontSize',12);
ylabel('Offset [km]','FontSize',12);

legend({'\Deltax','\Deltay','\Deltaz','Total'},'FontSize',10);


print(plot_core_cm, '-dpsc', 'CoreOffset.eps');















