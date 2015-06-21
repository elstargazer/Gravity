ccc

%% Load gravity model

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

Rref2=293000;

lmcosi_fa=plm2pot(lmcosi_grav,Rref2,mu,Rref,3,'nothing');

%% Load topography model

MaxDegreeTopo=20;
GravityFieldMaxDeg=18;
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% lmcosi_shape(4,3)=0;
[ri,~,~,~]=plm2xyz(lmcosi_shape,1,[],[],[]);

MeanRadius=lmcosi_shape(1,3);

ri=ri-MeanRadius;

%% Homogeneous gravity
Rref=265000;
mu=17.29e9;
rho=3457;
Resolution=1;
MaxDegreeTopo=100;
MaxDegreeGrav=20;
MaxTopoPower=8;

[ri2,~,~,~]=plm2xyz(lmcosi_shape,Resolution);

lmcosi_grav_shape=TopoSH2GravitySH(ri2,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_grav_shape(1,:)=[];

lmcosi_fa_shape=plm2pot(lmcosi_grav_shape,Rref2,mu,Rref,3,'nothing');

%% Computing free-air anomaly

[fa,~,~,~]=plm2xyz(lmcosi_fa,1,[],[],[]);
[fa_shape,~,~,~]=plm2xyz(lmcosi_fa_shape,1,[],[],[]);


%% Plotting Free-air anomaly

% [fa_fig,fa_ax]=MapRadialGrid(flipud(fa)*1e5);

MapRadialGrid(flipud(ri))
 
%% Global spectral power and admittance

Z_global=SphericalHarmonicAdmittance(lmcosi_fa,lmcosi_shape);

Z_global_shape=SphericalHarmonicAdmittance(lmcosi_fa_shape,lmcosi_shape);

Sgt_global=CrossPower(lmcosi_fa,lmcosi_shape);
Stt_global=CrossPower(lmcosi_shape,lmcosi_shape); 


Sgg_global=CrossPower(lmcosi_fa,lmcosi_fa);
Sgg_global_shape=CrossPower(lmcosi_fa_shape,lmcosi_fa_shape);
Sgg_global_cross_shape=CrossPower(lmcosi_fa,lmcosi_fa_shape);

% Z2_global=Sgt_global./Stt_global;



% plot(1:numel(Z2_global),Z2_global,'--or');

%% Global correlation

R_global=Sgg_global_cross_shape./sqrt(Sgg_global_shape.*Sgg_global);


%% Choosing region
% AGUaxes
% Boundary=inputm();
% 
% Boundary=fliplr(Boundary)*180/pi;

L=5;
L2=20;

% [V,C,jk1,jk2,XYZ]=localization(L,Boundary,[]);

fi=-90:5:90;
lambda=-180:5:180;

[fii,lambdai]=meshgrid(fi,lambda);

j=1;

circle_rad=50;
density=zeros(size(fii));

% progressbar('computing density')

% figure;
% 
% hold on;

progressbar(0);

GoodDegrees=2+L:GravityFieldMaxDeg-L;


for j=1:numel(density)

% patch coordinates
lat_center=fii(j);
lon_center=lambdai(j);

%  lat_center=-90;
%  lon_center=0;

colat_center=90-lat_center;


%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,lon_center,90-lat_center);


s=size(G2)


for i=1:s(2)
       
    lmcosi_window{i}=glm2lmcosi(G2,i);
    [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
end

%% Using localization

%[V,C,jk1,jk2,XYZ,~,G]=localization(L,'patch',[90-lat_center,lon_center,circle_rad]/180*pi);

%N=numel(C);

%lmcosi_fa_w=cell(size(C));
%lmcosi_shape_w=cell(size(C));

for i=1:s(2)
%    [r(:,:,i),lor,lar,Plm]=plm2xyz([jk1 jk2 C{i}],1);    
    
    
    [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,L2,'im');
    [lmcosi_fa_shape_w{i},~]=xyz2plm(r(:,:,i).*fa_shape,L2,'im');
    [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,L2,'im');
    
%     Z_local(:,i)=SphericalHarmonicAdmittance(lmcosi_fa_w{i},lmcosi_shape_w{i});
%     Z_local_shape(:,i)=SphericalHarmonicAdmittance(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
     
    Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
    Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    

    
    Sgt_local_shape(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
    Stt_local_shape(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i}); 
    
    
    Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
    Sgg_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_shape_w{i});
    Sgg_cross_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_w{i});    
   
end

% Z2_local=Sgt_local./Stt_local;

%% Localized admittance plot

Tapers=find(V2>0.99);
% Tapers=10:36;
% Tapers=1:5;

N0=(L+1)*(circle_rad/180);

% K=2;
% (K+1)*pi/(circle_rad/180*pi)-1

% plot_local=plot(Z_local(:,1:MaxTaper),'-.r');
% plot_local_shape=plot(Z_local_shape(:,1:MaxTaper),'--.r');

Sgt_local_mean=mean(Sgt_local(:,Tapers),2);
Stt_local_mean=mean(Stt_local(:,Tapers),2);

Sgt_local_mean_shape=mean(Sgt_local_shape(:,Tapers),2);
Stt_local_mean_shape=mean(Stt_local_shape(:,Tapers),2);

Z_local_mean=Sgt_local_mean./Stt_local_mean;
Z_local_mean_shape=Sgt_local_mean_shape./Stt_local_mean_shape;

Z_local_std=std(Sgt_local(:,Tapers)./Stt_local(:,Tapers),0,2);
Z_local_std_shape=std(Sgt_local_shape(:,Tapers)./Stt_local_shape(:,Tapers),0,2);

Degree=(1:numel(Z_local_mean))';

%% Density plot



density(j)=mean(Z_local_mean(GoodDegrees)./Z_local_mean_shape(GoodDegrees)*3457);

% density_std(j)=std(Z_local_mean(GoodDegrees)./Z_local_mean_shape(GoodDegrees)*3457);

% pcolor(density);
% 
% caxis([min(density(density~=0))-1 max(density(density~=0))+1]); 
% colorbar;

progressbar(j/numel(fii));
% 
% if ((density(j)<2000) || (density(j)>4000))
%     
%     disp([lambdai(j) fii(j)]); 
%     
% end

j/numel(fii)

end

AGUaxes
pcolorm(fii'/180*pi,lambdai'/180*pi,density2'); shading interp;





% contourm(fii/180*pi,lambdai/180*pi,density,20,'LineWidth',3); 

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



