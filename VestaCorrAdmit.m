ccc

FileName1='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
FileName2='/Users/antonermakov/Dawn/Gravity/VESTA20E/JGV20E02.SHA';
FileName3='/Users/antonermakov/Dawn/Gravity/VESTA10/JGV10C02.SHA';
FileName4='/Users/antonermakov/Dawn/Gravity/VESTA6/JGV06B02.SHA';

%% Load shape model
load VestaHASTALAVESTAshape_sh720.mat;
MaxTopoPower=7;
MaxDegreeTopo=300;
MaxDegreeGrav=20;

lmcosi=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);


[lmcosi_obs1,Rref,mu,mu_std]=ReadGRAILGravityModel(FileName1);
[lmcosi_obs2,Rref,mu,mu_std]=ReadGRAILGravityModel(FileName2);
[lmcosi_obs3,Rref,mu,mu_std]=ReadGRAILGravityModel(FileName3);
[lmcosi_obs4,Rref,mu,mu_std]=ReadGRAILGravityModel(FileName4);

%% Compute gravity from shape

ro=3.458942584995802e+03;
G=6.67384e-11;
mu=17.2883080293e9;
Rref=265000;


% making smooth Vesta
lmcosi_smooth=lmcosi;
lmcosi_smooth(2:end,3:4)=lmcosi_smooth(2:end,3:4)/5;


Resolution=0.25;

[ri_shape,~]=plm2xyz(lmcosi,Resolution);

lmcosi_gt=TopoSH2GravitySH(ri_shape,mu,ro,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_gt(1,:)=[];
lmcosi(1,:)=[];


[ri_shape_smooth,~]=plm2xyz(lmcosi_smooth,Resolution);

lmcosi_gt_smooth=TopoSH2GravitySH(ri_shape_smooth,mu,ro,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_smooth(1,:)=[];
lmcosi_gt_smooth(1,:)=[];




% SphericalHarmonicCorrelation(lmcosi,lmcosi_obs1,'m','--');
% SphericalHarmonicCorrelation(lmcosi,lmcosi_obs2,'b','--');
% SphericalHarmonicCorrelation(lmcosi,lmcosi_obs3,'g','--');
% SphericalHarmonicCorrelation(lmcosi,lmcosi_obs4,'r','--');
% 
% % legend({'HAMO2 (vesta20g)','LAMO (vesta20f)','HAMO (vesta10c)','Survey (vesta6b)'})
% 
% SphericalHarmonicCorrelation(lmcosi_gt,lmcosi_obs1,'m','-');
% SphericalHarmonicCorrelation(lmcosi_gt,lmcosi_obs2,'b','-');
% SphericalHarmonicCorrelation(lmcosi_gt,lmcosi_obs3,'g','-');
% SphericalHarmonicCorrelation(lmcosi_gt,lmcosi_obs4,'r','-');
% 
% legend({'HAMO2 (vesta20g)','LAMO (vesta20f)','HAMO (vesta10c)','Survey (vesta6b)'},'Location','SouthEast')
% 
% % SphericalHarmonicAdmittance(lmcosi_obs1,lmcosi,'c');
% % SphericalHarmonicAdmittance(lmcosi_obs2,lmcosi,'m');
% % SphericalHarmonicAdmittance(lmcosi_obs3,lmcosi,'y');
% % SphericalHarmonicAdmittance(lmcosi_obs4,lmcosi,'k');
% % 
% % legend({'JGV20F02','JGV20E02','JGV10C02','JGV06B02'})
% 
% 
% 
% set(gca,'XTick',2:2:20);
% xlabel('Degree','FontSize',12);
% ylabel('Correlation coefficient','FontSize',12);
% 
% axis([2 20 0 1]);
% 
% set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
% set(gcf, 'PaperPositionMode','auto')
% print(gcf, '-dpsc', 'GravityTopographyCorrelation.eps');

r=SphericalHarmonicCorrelation(lmcosi,lmcosi_gt,'m','-');
r_smooth=SphericalHarmonicCorrelation(lmcosi_smooth,lmcosi_gt_smooth,'m','-');



figure; hold on;
set(gca,'FontSize',20);

plot(r,'k-','LineWidth',3);
plot(r_smooth,'r-','LineWidth',3);

legend({'real shape','more spherical shape'});

axis([2 MaxDegreeGrav 0 1]);

xlabel('Degree','FontSize',20);
ylabel('Correlation','FontSize',20);




SurfRadialGrid(ri_shape)
axis equal
view(180+45,0);


SurfRadialGrid(ri_shape_smooth)
axis equal
view(180+45,0);




