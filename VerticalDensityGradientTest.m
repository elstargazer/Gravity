ccc

h=30000; % crustal thickness

dengrad=0:2:50;
densurf=2000:10:3000;

dengrad=dengrad/1000;

N=5;

G=6.67e-11;

MaxDegreeGrav=20;
MaxDegreeTopo=100;
MaxTopoPower=10;

rcore=110000;
fcore=0.1;
dencore=7800;
denmean=3457;

%% Doing something

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

[lmcosi_grav_obs,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

m=mu/G;

lmcosi_grav_obs=AddZeroHarm(lmcosi_grav_obs,1);

acore=rcore/((1-fcore).^(1/3));
ccore=rcore/((1-fcore).^(1/3))-fcore*rcore/((1-fcore).^(1/3));

Ccore=SHRotationalEllipsoid(acore,ccore,MaxDegreeGrav,MaxDegreeGrav,Rref);
lmcosi_core=CS2lmcosi(Ccore,Ccore*0);
lmcosi_core=AddZeroHarm(lmcosi_core,1);

Vcore=4/3*pi*rcore^3;

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

[dengradi,densurfi]=meshgrid(dengrad,densurf);

% figure; hold on;
n=1:MaxDegreeGrav;

GoodDegrees=3:17;


for j=1:numel(dengradi)
    
    

    den=densurfi(j)+dengradi(j).*z;

    lmcosi_grav_total=lmcosi_grav_homo;
    lmcosi_grav_total(:,3:4)=0;

    V=4/3*pi*R.^3;

    Vmantle=V(end);

    diffden=diff(den(1:end-1));

    denlayer=[densurfi(j) diffden];

    mlayer=V(1:end-1).*denlayer;

    mtotal=sum(mlayer);

    Vinner=V(end);

    denmantlelayer=-((m-mtotal)-dencore*Vcore+sum(denlayer)*Vcore)/(Vcore-Vmantle);


    if (denmantlelayer<0)
        d2(j)=NaN;
    else 

    

    dencorelayer=((m-mtotal)-dencore*Vmantle+sum(denlayer)*Vmantle)/(Vcore-Vmantle);

    mmantle=denmantlelayer*Vmantle;
    mcore=dencorelayer*Vcore;



    for i=1:numel(z)-1
    
        lmcosi_grav_total(:,3:4)=lmcosi_grav_total(:,3:4)+denlayer(i)*lmcosi_grav{i}(:,3:4)/mu;
    
    end


    %% Adding mantle and the core

    lmcosi_grav_total(:,3:4)=lmcosi_grav_total(:,3:4)+denmantlelayer*lmcosi_grav{end}(:,3:4)/mu;

    lmcosi_grav_total(:,3:4)=lmcosi_grav_total(:,3:4)+lmcosi_core(:,3:4)*mcore/m;

    denstructurelayer=[denlayer denmantlelayer dencorelayer];

    cumsum(denstructurelayer);

    %% Plotting density structure



    %% Computing admittance

    Adm=SphericalHarmonicAdmittance(lmcosi_grav_total,lmcosi_shape);
    Adm_homo=SphericalHarmonicAdmittance(lmcosi_grav_homo,lmcosi_shape);
    Adm_obs=SphericalHarmonicAdmittance(lmcosi_grav_obs,lmcosi_shape);

    %% Computing effective density

    deneff=Adm./Adm_homo.*denmean;
    deneff_obs=Adm_obs./Adm_homo.*denmean;


    %% chi squared

    d2(j)=sum((deneff(GoodDegrees)-deneff_obs(GoodDegrees)).^2);


% plot(n,deneff,'-b','LineWidth',3);
% plot(n,deneff_obs,'-r','LineWidth',3);
% drawnow;
% pause(0.05);

% unplot(2);
    end

end

d2=reshape(d2,size(densurfi));

figure; hold on;
set(gca,'FontSize',20);


pcolor(densurfi,dengradi*1000,log10(d2)); shading interp;
contour(densurfi,dengradi*1000,log10(d2),50,'Color','k');

index=find(d2==min(d2(:)));

plot(densurfi(index),dengradi(index)*1000,'w*','MarkerSize',20);

cbar=colorbar;

xlabel('Surface density [kg/m^{3}]','FontSize',20);
ylabel('Vertical density gradient [kg/m ^{3}/km]','FontSize',20);

ylabel(cbar,'log_{10}(\Delta^{2})','FontSize',20);

box on;


%% Plotting effective density spectrum


n=1:MaxDegreeGrav;

figure; hold on;
set(gca,'FontSize',20);

plot(n,deneff,'-b','LineWidth',3);
plot(n,deneff_obs,'-r','LineWidth',3);

xlabel('Degree','FontSize',20);
ylabel('Effective density [kg/m^{3}]','FontSize',20);

box on;


%% Plot Density structure

den=densurfi(index)+dengradi(index).*z;
V=4/3*pi*R.^3;
Vmantle=V(end);
diffden=diff(den(1:end-1));
denlayer=[densurfi(index) diffden];
mlayer=V(1:end-1).*denlayer;
mtotal=sum(mlayer);
Vinner=V(end);
denmantlelayer=-((m-mtotal)-dencore*Vcore+sum(denlayer)*Vcore)/(Vcore-Vmantle);
dencorelayer=((m-mtotal)-dencore*Vmantle+sum(denlayer)*Vmantle)/(Vcore-Vmantle);
mmantle=denmantlelayer*Vmantle;
mcore=dencorelayer*Vcore;
denstructurelayer=[denlayer denmantlelayer dencorelayer];
denstructurecum=cumsum(denstructurelayer);

figure; hold on;
set(gca,'FontSize',20);

for i=1:numel(z)-1    
    
    plot(R(1)-linspace(z(i),z(i+1),10),linspace(denstructurecum(i),denstructurecum(i),10),'-k','LineWidth',3) 
    plot(R(1)-linspace(z(i+1),z(i+1),10),linspace(denstructurecum(i),denstructurecum(i+1),10),'-k','LineWidth',3)
    
end

plot(linspace(rcore,R(1)-z(end),10),linspace(denstructurecum(i+1),denstructurecum(i+1),10),'-k','LineWidth',3);
plot(linspace(rcore,rcore,10),linspace(denstructurecum(i+1),dencore,10),'-k','LineWidth',3)
plot(linspace(0,rcore,10),linspace(dencore,dencore,10),'-k','LineWidth',3);

xlim([0 R(1)])

xlabel('Mean radius [m]','FontSize',20);
ylabel('Density [kg/m^{3}]','FontSize',20);

box on






















