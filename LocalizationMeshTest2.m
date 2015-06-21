ccc
%% Load gravity model

GravityFileName='~/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);


Rref2=293000;

lmcosi_fa=plm2pot(lmcosi_grav,Rref2,mu,Rref,3,'nothing');
% [lmcosi_fa,~,~]=plm2rot(lmcosi_fa,0,90,0,'dlmb');
%% Load topography model

MaxDegreeTopo=20;
GravityFieldMaxDeg=18;
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% [lmcosi_shape,~,~]=plm2rot(lmcosi_shape,0,90,0,'dlmb');

% lmcosi_shape(4,3)=0;
[ri,~,~,~]=plm2xyz(lmcosi_shape,1,[],[],[]);

MeanRadius=lmcosi_shape(1,3);

ri=ri-MeanRadius;

MapRadialGrid(flipud(ri));

%% Homogeneous gravity
Rref=265000;
mu=17.29e9;
rho=3457;
rho_mean=rho;
Resolution=1;
MaxDegreeTopo=100;
MaxDegreeGrav=20;
MaxTopoPower=10;

[ri2,~,~,~]=plm2xyz(lmcosi_shape,Resolution);

lmcosi_grav_shape=TopoSH2GravitySH(ri2,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_grav_shape(1,:)=[];

lmcosi_fa_shape=plm2pot(lmcosi_grav_shape,Rref2,mu,Rref,3,'nothing');

%% Computing free-air anomaly

[fa,~,~,~]=plm2xyz(lmcosi_fa,1,[],[],[]);
[fa_shape,~,~,~]=plm2xyz(lmcosi_fa_shape,1,[],[],[]);


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

%% Global correlation

R_global=Sgg_global_cross_shape./sqrt(Sgg_global_shape.*Sgg_global);


%% Choosing geo mesh

L=5;
L2=20;

fi=-90:5:90;
lambda=-180:5:180;

[fii,lambdai]=meshgrid(fi,lambda);


circle_rad=50;
density_mean=zeros(size(fii));
p1=zeros(size(fii));
p2=zeros(size(fii));


progressbar(0);

GoodDegrees=2+L:GravityFieldMaxDeg-L;

%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,0,0);

s=size(G2);


for i=1:s(2)       
    lmcosi_window_basic{i}=glm2lmcosi(G2,i);    
end

Tapers=find(V2>0.75);
NTapers=numel(Tapers);
lmcosi_window_basic=lmcosi_window_basic(Tapers);


for j=1:numel(density_mean)

% patch coordinates
lat_center=fii(j);
lon_center=lambdai(j);

colat_center=90-lat_center;


alp=0;
bta=lat_center-90; 
gam=-lon_center;

for i=1:NTapers
    
    [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,bta,gam,'dlmb');
    
    [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
    
    
    [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,L2,'im');
    [lmcosi_fa_shape_w{i},~]=xyz2plm(r(:,:,i).*fa_shape,L2,'im');
    [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,L2,'im');
         
    Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
    Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    

    
    Sgt_local_shape(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
    Stt_local_shape(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i}); 
    
    
    Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
    Sgg_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_shape_w{i});
    Sgg_cross_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_w{i});    
   
end


%% Localized admittance plot

N0=(L+1)*(circle_rad/180);

Sgt_local_mean=mean(Sgt_local,2);
Stt_local_mean=mean(Stt_local,2);

Sgt_local_mean_shape=mean(Sgt_local_shape,2);
Stt_local_mean_shape=mean(Stt_local_shape,2);

Z_local_mean=Sgt_local_mean./Stt_local_mean;
Z_local_mean_shape=Sgt_local_mean_shape./Stt_local_mean_shape;

Z_local_std=std(Sgt_local./Stt_local,0,2);
Z_local_std_shape=std(Sgt_local_shape./Stt_local_shape,0,2);

Degree=(1:numel(Z_local_mean))';


Z_local_mean_std=Z_local_std/sqrt(NTapers);
Z_local_mean_std_shape=Z_local_std_shape/sqrt(NTapers);



drhodx=rho_mean./Z_local_mean_shape;
drhody=-rho_mean*Z_local_mean./(Z_local_mean_shape.^2);

density_local_std=sqrt((drhodx.^2).*(Z_local_mean_std.^2)+(drhody.^2).*(Z_local_mean_std_shape.^2))

%% Density 

density=Z_local_mean(GoodDegrees)./Z_local_mean_shape(GoodDegrees)*rho;

p=polyfit(GoodDegrees,density',1);

p1(j)=p(1);
p2(j)=p(2);

density_mean(j)=mean(density);



progressbar(j/numel(fii));

j/numel(fii)

end

density_mean=reshape(density_mean,size(fii));

AGUaxes
pcolorm(fii'/180*pi,lambdai'/180*pi,density_mean'); shading interp;
colorbar
PlotBasins

AGUaxes
pcolorm(fii'/180*pi,lambdai'/180*pi,p1'); shading interp;
colorbar
PlotBasins

AGUaxes
pcolorm(fii'/180*pi,lambdai'/180*pi,p2'); shading interp;
colorbar
PlotBasins

% AGUaxes
% pcolorm(fii'/180*pi,lambdai'/180*pi,abs(p2./p1)'); shading interp;
% colorbar


density_middle=zeros(size(density_mean));

for i=1:numel(fii)
    
p=[p1(i) p2(i)]; 
density_middle(i)=polyval(p,mean(GoodDegrees));
    
end
% 
% AGUaxes
% pcolorm(fii'/180*pi,lambdai'/180*pi,density_middle); shading interp;
% colorbar

figure
plot(p1(:),density_middle(:),'.')
xlabel('Slope','FontSize',20);
ylabel('Intercept','FontSize',20);
set(gca,'FontSize',20)
xlabel('Slope [kg/m^3/degree]','FontSize',20);
ylabel('Intercept [kg/m^3]','FontSize',20);




%% Plot Rheasilvia and Veneniea

R_R=450/2/265*180/pi;
R_V=395/2/265*180/pi;

fi_R=-75;
lambda_R=301;
fi_V=-52;
lambda_V=170;

[latc_R,lonc_R] = scircle1(fi_R,lambda_R,R_R);
[latc_V,lonc_V] = scircle1(fi_V,lambda_V,R_V);

plotm(latc_R/180*pi,lonc_R/180*pi,'--k','LineWidth',4);
plotm(latc_V/180*pi,lonc_V/180*pi,'--k','LineWidth',4);



