ccc

%% Input Paramerers;

% GravityFileName='/Users/antonermakov/GRAIL/Gravity_Models/jggrx_0660b_sha.tab.txt';
GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

MaxDegreeTopo=100;
MaxDegreeGrav=20;
Rref_gravity=293000;
% Rref_gravity=1738000;
rho=3457;
Resolution=1;
MaxTopoPower=7;
L=5;
MinConcentration=0.70;
NTess=2;
circle_rad=50;
MaxDegreeExp=10;
% rho_mean=3457.5;
GoodDegrees=3+L:MaxDegreeGrav-L-2;

%% Load gravity model

[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

% load lmcosi_calc_new.mat
% lmcosi_grav=lmcosi_calc_new;
% lmcosi_grav=TruncateGravityModel(lmcosi_grav,MaxDegreeGrav,0);

% lmcosi_grav(3,3)=0;

lmcosi_fa=plm2pot(lmcosi_grav,Rref_gravity,mu,Rref,3,'nothing');

%% Load topography model

load VestaHASTALAVESTAshape_sh720.mat
% lmcosi_shape=load('/Users/antonermakov/GRAIL/Topography/MoonTopo2600p.shape.sh');

lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% lmcosi_grav(4,3)=0;

[ri,~,~,~]=plm2xyz(lmcosi_shape,Resolution);
MeanRadius=lmcosi_shape(1,3);
ri=ri-MeanRadius;
MapRadialGrid(flipud(ri));

PlotOlivine

%% Homogeneous gravity

lmcosi_grav_shape=TopoSH2GravitySH(ri+MeanRadius,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);
% [U,~,~,~]=plm2xyz(lmcosi_grav_shape,Resolution,[],[],[]);

lmcosi_grav_shape(1,:)=[];
lmcosi_fa_shape=plm2pot(lmcosi_grav_shape,Rref_gravity,mu,Rref,3,'nothing');

%% Computing free-air anomaly

[fa,~,~,~]=plm2xyz(lmcosi_fa,Resolution,[],[],[]);
[fa_shape,~,~,~]=plm2xyz(lmcosi_fa_shape,Resolution,[],[],[]);


%% Plotting Free-air anomaly

% [fa_fig,fa_ax]=MapRadialGrid(flipud(fa)*1e5);
% MapRadialGrid(flipud(ri))
 
%% Global spectral power and admittance

Z_global=SphericalHarmonicAdmittance(lmcosi_fa,lmcosi_shape);
Z_global_shape=SphericalHarmonicAdmittance(lmcosi_fa_shape,lmcosi_shape);

Sgt_global=CrossPower(lmcosi_fa,lmcosi_shape);
Stt_global=CrossPower(lmcosi_shape,lmcosi_shape); 

Sgg_global=CrossPower(lmcosi_fa,lmcosi_fa);
Sgg_global_shape=CrossPower(lmcosi_fa_shape,lmcosi_fa_shape);
Sgg_global_cross_shape=CrossPower(lmcosi_fa,lmcosi_fa_shape);

figure; hold on;
set(gca,'FontSize',20);

degree=1:MaxDegreeGrav;
plot(degree,Z_global,'-ok','LineWidth',4,'MarkerSize',4);
plot(degree,Z_global_shape,'-or','LineWidth',4,'MarkerSize',4);

%% Global correlation
R_global=Sgg_global_cross_shape./sqrt(Sgg_global_shape.*Sgg_global);

%% Icosahedron mesh;
TR=IcosahedronMesh;
TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;
x_t=TR_2.X(:,1);
y_t=TR_2.X(:,2);
z_t=TR_2.X(:,3);

[lambdai,fii,~]=cart2sph(x_t,y_t,z_t);
lambdai=lambdai*180/pi;
fii=fii*180/pi;

density_mean=zeros(size(fii));
p1=zeros(size(fii));
p2=zeros(size(fii));

Npoints=numel(fii);
% MaxDegreeExp=fix(0.5*(-3+sqrt(1+8*Npoints)))

%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,0,0);

s=size(G2);

for i=1:s(2)       
    lmcosi_window_basic{i}=glm2lmcosi(G2,i);    
end

Tapers=find(V2>MinConcentration);
NumberOfTapers=numel(Tapers);
lmcosi_window_basic=lmcosi_window_basic(Tapers);

NTapers=sum(V2>MinConcentration);

progressbar(0);

for j=1:numel(density_mean)

% patch coordinates
lat_center=fii(j);
lon_center=lambdai(j);

colat_center=90-lat_center;

alp=0;
bta=lat_center-90; 
gam=-lon_center;

for i=1:NumberOfTapers
    
    [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,bta,gam,'dlmb');
    
    [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},Resolution);
    
    
    [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,MaxDegreeGrav,'im');
    [lmcosi_fa_shape_w{i},~]=xyz2plm(r(:,:,i).*fa_shape,MaxDegreeGrav,'im');
    [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,MaxDegreeGrav,'im');
         
    Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
    Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    

    
    Sgt_local_shape(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
    Stt_local_shape(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i}); 
        
    Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
    Sgg_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_shape_w{i});
    Sgg_cross_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_w{i});    
    
    Corr(:,i)=SphericalHarmonicCorrelation(lmcosi_fa_shape_w{i},lmcosi_fa_w{i});
    
    
%     q_topo(:,i)=IsotropicRatio(lmcosi_Moon_shape(2:end,:),lmcosi_Moon_shape(2:end,:));
    q_grav(:,i)=IsotropicRatio(lmcosi_fa_w{i},lmcosi_fa_w{i});
   
end

%% Localized admittance plot

N0=(L+1)*(circle_rad/180);


q_grav_mean(:,j)=mean(q_grav,2);
q_grav_std(:,j)=std(q_grav,0,2);


Sgt_local_mean=mean(Sgt_local,2);
Stt_local_mean=mean(Stt_local,2);

Sgt_local_mean_shape=mean(Sgt_local_shape,2);
Stt_local_mean_shape=mean(Stt_local_shape,2);

Z_local_mean=Sgt_local_mean./Stt_local_mean;
Z_local_mean_shape=Sgt_local_mean_shape./Stt_local_mean_shape;

Z_local_std=std(Sgt_local./Stt_local,0,2);
Z_local_std_shape=std(Sgt_local_shape./Stt_local_shape,0,2);

Z_local_mean_std=Z_local_std/sqrt(NTapers);
Z_local_mean_std_shape=Z_local_std_shape/sqrt(NTapers);

drhodx=rho./Z_local_mean_shape;
drhody=-rho*Z_local_mean./(Z_local_mean_shape.^2);

density_local_std=(1/sqrt(numel(GoodDegrees)))*sqrt((drhodx.^2).*(Z_local_mean_std.^2)+(drhody.^2).*(Z_local_mean_std_shape.^2));

Degree=(1:numel(Z_local_mean))';

Z_t(:,j)=Z_local_mean(GoodDegrees);
Z_t_shape(:,j)=Z_local_mean_shape(GoodDegrees);

Corr_mean(:,j)=mean(Corr,2);

%% Density 
density=Z_local_mean(GoodDegrees)./Z_local_mean_shape(GoodDegrees)*rho;
density_t(:,j)=density;

p=polyfit(GoodDegrees,density',1);

[a,b,~,~,~,~,~]=wtls_line(GoodDegrees,density',...
    ones(1,numel(GoodDegrees)),density_local_std(GoodDegrees)');

density_wlf=a*GoodDegrees+b;

% errorbar(GoodDegrees,density',density_local_std(GoodDegrees));

p1(j)=p(1);
p2(j)=p(2);


density_mean(j)=mean(density_wlf);
% density_mean(j)=mean(density);


density_mean_w(j)=sum(density.*(1./density_local_std(GoodDegrees).^2))./sum(1./density_local_std(GoodDegrees).^2);

% fitobject = fit(GoodDegrees,y,fitType,...,'Weight', Weights)

progressbar(j/numel(fii));

j/numel(fii)

end



density_mean=reshape(density_mean,size(fii));
density_mean_w=reshape(density_mean_w,size(fii));

% density_mean=density_mean_w;



lmcosi_den=xyz2plm(density_mean,MaxDegreeExp,'irr',fii,lambdai);
[density_mean_sh,lonsh,latsh]=plm2xyz(lmcosi_den,.2);
[lonsh,latsh]=meshgrid(lonsh,latsh);

AGUaxes
pcolorm(latsh/180*pi,lonsh/180*pi,density_mean_sh); shading interp
cb1=colorbar('FontSize',20);
ylabel(cb1,'Effective density [kg/m ^3] ','FontSize',20);
PlotVestaFeatures
title('Effective density ','FontSize',20);

PlotOlivine


WriteXYZ(lonsh,latsh,density_mean_sh,'MeanEffectiveDensity.xyz');


lmcosi_p1=xyz2plm(p1,MaxDegreeExp,'irr',fii,lambdai);
[p1_sh,lonsh,latsh]=plm2xyz(lmcosi_p1,.2);
[lonsh,latsh]=meshgrid(lonsh,latsh);

AGUaxes
pcolorm(latsh/180*pi,lonsh/180*pi,p1_sh); shading interp
cb1=colorbar('FontSize',20);
ylabel(cb1,'Effective density slope [kg/m ^3/degree] ','FontSize',20);
PlotVestaFeatures
title('Effective density slope ','FontSize',20);

WriteXYZ(lonsh,latsh,p1_sh,'EffectiveDensitySlope.xyz');


AGUaxes
scatterm(fii/180*pi,lambdai/180*pi,400,density_mean,'filled'); 
cb1=colorbar('FontSize',20);
ylabel(cb1,'Effective density [kg/m ^3] ','FontSize',20);
PlotVestaFeatures
title('Effective density ','FontSize',20);

for i=1:numel(fii)    
    [latc,lonc] = scircle1(fii(i),lambdai(i),circle_rad);
    plotm(latc/180*pi,lonc/180*pi,'-r','LineWidth',4);       
end

AGUaxes
scatterm(fii/180*pi,lambdai/180*pi,200,p1,'filled'); 
cb2=colorbar('FontSize',20);
ylabel(cb2,'Effective density slope [kg/m ^3/degree] ','FontSize',20);
PlotVestaFeatures
title('Effective density slope ','FontSize',20);


WriteXYZ(lambdai,fii,density_mean,'IsoMeshAdmittance.xyz');

% AGUaxes
% pcolorm(fii'/180*pi,lambdai'/180*pi,abs(p2./p1)'); shading interp;
% colorbar

% density_middle=zeros(size(density_mean));
% 
% for i=1:numel(fii)
%     
% p=[p1(i) p2(i)]; 
% density_middle(i)=polyval(p,mean(GoodDegrees));
%     
% end
% 
% AGUaxes
% pcolorm(fii'/180*pi,lambdai'/180*pi,density_middle); shading interp;
% colorbar

Npoints=100000;
[fi_rand,lambda_rand]=GenerateRandomSphericalCoord(Npoints);

p1_rand=griddata(fii,lambdai,p1,fi_rand*180/pi,lambda_rand*180/pi,'nearest');
density_mean_rand=griddata(fii,lambdai,density_mean,fi_rand*180/pi,lambda_rand*180/pi,'nearest');

figure; hold on;
plot(p1(:),density_mean(:),'.','MarkerSize',18)
set(gca,'FontSize',20)
xlabel('Density variation per SH degree [kg/m^3/degree] ','FontSize',20);
ylabel('Effective density [kg/m ^3] ','FontSize',20);

pa=polyfit(p1(:),density_mean(:),1);
s=min(p1(:)):1:max(p1(:));
densitya=polyval(pa,s);
plot(s,densitya,'-r','LineWidth',3);

figure; hold on;
plot(p1_rand(:),density_mean_rand(:),'.','MarkerSize',18)
set(gca,'FontSize',20)
xlabel('Density variation per SH degree [kg/m^3/degree]','FontSize',20);
ylabel('Effective density [kg/m ^3] ','FontSize',20);

pa=polyfit(p1_rand(:),density_mean_rand(:),1);
s=min(p1_rand(:)):1:max(p1_rand(:));
densitya=polyval(pa,s);
plot(s,densitya,'-r','LineWidth',3);


fitobject = fit(p1_rand(:),density_mean_rand(:),'poly1');

fit_pl=plot(fitobject,'m')

set(fit_pl,'LineWidth',3);


xlabel('Density variation per SH degree [kg/m^3/degree]','FontSize',20);
ylabel('Effective density [kg/m ^3] ','FontSize',20);



%% Isoptropic ratio


% lmcosi_q=cell(1,numel(GoodDegrees));
AGUaxes
PlotVestaFeatures
caxis([0.8 4]);


vidObj = VideoWriter('IsotropicRatioMovie.avi');
open(vidObj);

for i=GoodDegrees(1):GoodDegrees(end)

% for i=1:numel(GoodDegrees)
    
    lmcosi_q=xyz2plm(q_grav_mean(i,:),MaxDegreeExp,'irr',fii,lambdai);
    [q_grav_sh,lonsh,latsh]=plm2xyz(lmcosi_q,.2);
    [lonsh,latsh]=meshgrid(lonsh,latsh);
    
% end

pcolorm(latsh/180*pi,lonsh/180*pi,q_grav_sh); shading interp
cb1=colorbar('FontSize',20);
ylabel(cb1,'Isotropic ratio []','FontSize',20);
title(['Isotropic ratio at degree ' num2str(i) ' '],'FontSize',20);


   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);

end


close(vidObj);

%% Clusters

NClucters=10;

[idx,ctrs] = kmeans(density_t',NClucters);
cc=jet(NClucters);


AGUaxes;
scatterm(fii/180*pi,lambdai/180*pi,100,idx,'filled')
PlotVestaFeatures;


% Voronoi diagram

xu=cos(fii/180*pi).*cos(lambdai/180*pi);
yu=cos(fii/180*pi).*sin(lambdai/180*pi);
zu=sin(fii/180*pi);

xyzu=[xu yu zu]';

[P, K, voronoiboundary] = voronoisphere(xyzu,'resolution', 0.1/180*pi);

% Graphic

f = figure(1);
clf(f);
set(f,'Renderer','zbuffer');
ax = axes('Parent', f);
hold(ax, 'on');
axis(ax,'equal');

plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');
clmap = cool();
ncl = size(clmap,1);


for k = 1:numel(xu)
    X = voronoiboundary{k};
    cl = clmap(mod(k,ncl)+1,:);
    fill3(X(1,:),X(2,:),X(3,:),cc(idx(k),:),'Parent',ax,'EdgeColor','w');
end
axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);

% fii,lambdai


% Mapping clusters
AGUaxes; hold on


% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');

ncl = size(clmap,1);

for k = 1:numel(xu)
    X = voronoiboundary{k};
    cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    fillm(fi_b,lambda_b,0,cc(idx(k),:),'EdgeColor','none');
end

PlotVestaFeatures

