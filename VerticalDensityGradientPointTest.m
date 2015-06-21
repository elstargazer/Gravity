ccc

h=30000;

dengrad=150;
densurf=1500;

dengrad=dengrad/1000;

N=5;

G=6.67e-11;

MaxDegreeGrav=20;
MaxDegreeTopo=100;
MaxTopoPower=10;


denmean=3457;

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

[lmcosi_grav_obs,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

lmcosi_grav_obs=AddZeroHarm(lmcosi_grav_obs,1);

rcore=110000;
fcore=0.1;

acore=rcore/((1-fcore).^(1/3));
ccore=rcore/((1-fcore).^(1/3))-fcore*rcore/((1-fcore).^(1/3));

Ccore=SHRotationalEllipsoid(acore,ccore,MaxDegreeGrav,MaxDegreeGrav,Rref);



z=0:h/N:h;


%% Load topograpy model


Resolution=1;

% int_str_fig=figure('color',[0 0 0]);
% plotplm(lmcosi_mantle_shape_new,[],[],2,1,[],[],[]);
% alpha(0.5)
% hold on
%  
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

[ri_shape,lambda,fi,~]=plm2xyz(lmcosi_shape,Resolution, [0 90 360 -90]);

[lambdai,fii]=meshgrid(lambda,fi);

[xi,yi,zi]=sph2cart(lambdai/180*pi,fii/180*pi,ri_shape);

% figure; hold on;
% 
% surf(xi(:,1:180),yi(:,1:180),zi(:,1:180),ri_shape(:,1:180));
% shading interp;
% alpha(0.5);

%% Homogeneous model

lmcosi_grav_homo=TopoSH2GravitySH(ri_shape,mu,denmean,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

%% Making crust


for i=1:numel(z)
    
    ri_shape_inner(:,:,i)=ri_shape-z(i);
    [xi,yi,zi]=sph2cart(lambdai/180*pi,fii/180*pi,ri_shape_inner(:,:,i));
    
%     surf(xi(:,1:180),yi(:,1:180),zi(:,1:180),ri_shape_inner(:,1:180,i));  shading interp; alpha(0.5);  
    
    lmcosi_shape_inner{i}=xyz2plm(ri_shape_inner(:,:,i),MaxDegreeTopo,'im');
    
    R(i)=lmcosi_shape_inner{i}(1,3);
    
end

% view(0,0);


%% Computing overall gravity coefficients

for i=1:numel(z)

lmcosi_grav{i}=TopoSH2GravitySH(ri_shape_inner(:,:,i),1,1,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

end


%% Main loop



den=densurf+dengrad.*z;

lmcosi_grav_total=lmcosi_grav_homo;
lmcosi_grav_total(:,3:4)=0;

V=4/3*pi*R.^3;

diffden=diff(den);

denlayer=[densurf diffden];

mulayer=G*V.*denlayer;

mutotal=sum(mulayer);

Vinner=V(end);

deninner=(mu-mutotal)/G/Vinner;

denlayer(end)=deninner;

mulayer(end)=mulayer(end)+(mu-mutotal);




for i=1:numel(z)
    
    lmcosi_grav_total(:,3:4)=lmcosi_grav_total(:,3:4)+denlayer(i)*lmcosi_grav{i}(:,3:4)/mu;
    
end


%% Adding mantle and the core



%% Plotting density structure



%% Computing admittance

Adm=SphericalHarmonicAdmittance(lmcosi_grav_total,lmcosi_shape);
Adm_homo=SphericalHarmonicAdmittance(lmcosi_grav_homo,lmcosi_shape);
Adm_obs=SphericalHarmonicAdmittance(lmcosi_grav_obs,lmcosi_shape);

%% Computing effective density

deneff=Adm./Adm_homo.*denmean;
deneff_obs=Adm_obs./Adm_homo.*denmean;


%% chi squared

GoodDegrees=3:17;


d2=sum((deneff(GoodDegrees)-deneff_obs(GoodDegrees)).^2);




%% Plotting effective density spectrum


n=1:MaxDegreeGrav;

figure; hold on;
set(gca,'FontSize',20);

plot(n,deneff,'-b','LineWidth',3);
plot(n,deneff_obs,'-r','LineWidth',3);

xlabel('Degree','FontSize',20);
ylabel('Effective density [kg/m^{3}]','FontSize',20);





