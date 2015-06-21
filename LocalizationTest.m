ccc
%
%% Load gravity model

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

Rref2=293000;

lmcosi_fa=plm2pot(lmcosi_grav,Rref2,mu,Rref,3,'nothing');

%% Load topography model

MaxDegreeTopo=20;

load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

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
MaxTopoPower=15;
rho_mean=3457;
MaxLambda=0.7;
circle_rad=60;

[ri2,~,~,~]=plm2xyz(lmcosi_shape,Resolution);


lmcosi_grav_shape=TopoSH2GravitySH(ri2,mu,rho,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_grav_shape(1,:)=[];

lmcosi_fa_shape=plm2pot(lmcosi_grav_shape,Rref2,mu,Rref,3,'nothing');

%% Computing free-air anomaly

[fa,~,~,~]=plm2xyz(lmcosi_fa,1,[],[],[]);
[fa_shape,~,~,~]=plm2xyz(lmcosi_fa_shape,1,[],[],[]);


%% Plotting Free-air anomaly

[fa_fig,fa_ax]=MapRadialGrid(flipud(fa_shape)*1e5);
 
%% Global spectral power and admittance

Z_global=SphericalHarmonicAdmittance(lmcosi_fa,lmcosi_shape);

Z_global_shape=SphericalHarmonicAdmittance(lmcosi_fa_shape,lmcosi_shape);

Sgt_global=CrossPower(lmcosi_fa,lmcosi_shape);
Stt_global=CrossPower(lmcosi_shape,lmcosi_shape); 


Sgg_global=CrossPower(lmcosi_fa,lmcosi_fa);
Sgg_global_shape=CrossPower(lmcosi_fa_shape,lmcosi_fa_shape);
Sgg_global_cross_shape=CrossPower(lmcosi_fa,lmcosi_fa_shape);

Adm_fig=figure; hold on;

plot_global=plot(1:numel(Z_global),1e5*1000*Z_global,'-ob','LineWidth',5);
plot_global_shape=plot(1:numel(Z_global),1e5*1000*Z_global_shape,'--ob','LineWidth',5);

xlim([1 numel(Z_global)]);

xlabel('Degree','FontSize',20);
ylabel('Admittance [mGal/km]','FontSize',20);

Z2_global=Sgt_global./Stt_global;

ylim(1e5*1000*[min(Z_global) max(Z_global)])

% plot(1:numel(Z2_global),Z2_global,'--or');

%% Global correlation

R_global=Sgg_global_cross_shape./sqrt(Sgg_global_shape.*Sgg_global);

Corr_fig=figure; hold on;

plot_global_corr=plot(1:numel(R_global),R_global,'-ob','LineWidth',5);

xlabel('Degree','FontSize',20);
ylabel('Correlation','FontSize',20);

set(gca,'FontSize',20);

%% Choosing region
% AGUaxes
% Boundary=inputm();
% 
% Boundary=fliplr(Boundary)*180/pi;

L=5;
L2=20;


% [V,C,jk1,jk2,XYZ]=localization(L,Boundary,[]);


% patch coordinates
lon_center=301;
lat_center=-75;


colat_center=90-lat_center;


[latc,lonc] = scircle1(lat_center,lon_center,circle_rad);

axes(fa_ax)
plotm(latc/180*pi,lonc/180*pi,'--k','LineWidth',4);

%% Using gplalpha
% [G2,V2,EL,EM,N2,GM2AL,MTAP,IMTAP]=glmalpha(circle_rad,L,1,0,0,0);
% 
% 
% 
% G2=rotateGp(G2,-40,145);
% 
% s=size(G2)
% 
% progressbar('computing windows');
% 
% for i=1:s(2)
%     
%     progressbar(i/s(2));
%     
%     lmcosi_window{i}=glm2lmcosi(G2,i);
%     [r2(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
% end



progressbar(1);

%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,lon_center,90-lat_center);

s=size(G2)

progressbar('computing windows');

for i=1:s(2)
    
    progressbar(i/s(2));
    
    lmcosi_window{i}=glm2lmcosi(G2,i);
    [r2(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
end



%% Using localization

%[V,C,jk1,jk2,XYZ,~,G]=localization(L,'patch',[90-lat_center,lon_center,circle_rad]/180*pi);

%N=numel(C);

%lmcosi_fa_w=cell(size(C));
%lmcosi_shape_w=cell(size(C));

r=r2;



for i=1:s(2)
%    [r(:,:,i),lor,lar,Plm]=plm2xyz([jk1 jk2 C{i}],1);    
    
    
    [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,L2,'im');
    [lmcosi_fa_shape_w{i},~]=xyz2plm(r(:,:,i).*fa_shape,L2,'im');
    [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,L2,'im');
    
    Z_local(:,i)=SphericalHarmonicAdmittance(lmcosi_fa_w{i},lmcosi_shape_w{i});
    Z_local_shape(:,i)=SphericalHarmonicAdmittance(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
    
    Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
    Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    

    
    Sgt_local_shape(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_shape_w{i});
    Stt_local_shape(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i}); 
    
    
    Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
    Sgg_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_shape_w{i});
    Sgg_cross_shape_local(:,i)=CrossPower(lmcosi_fa_shape_w{i},lmcosi_fa_w{i});
    
   
end

Z2_local=Sgt_local./Stt_local;


%% Eigenvalue plot = concentration factor
Eigenvalue_fig=figure; hold on;

plot(V2,'-.','LineWidth',3);

grid on

%% Localized admittance plot

figure(Adm_fig);

Tapers=find(V2>MaxLambda);

NTapers=sum(V2>MaxLambda);
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

plot_local_mean=plot(1e5*1000*Z_local_mean,'-or','LineWidth',4);
plot_local_mean_shape=plot(1e5*1000*Z_local_mean_shape,'--or','LineWidth',4);

Z_local_std=std(Sgt_local(:,Tapers)./Stt_local(:,Tapers),0,2);
Z_local_std_shape=std(Sgt_local_shape(:,Tapers)./Stt_local_shape(:,Tapers),0,2);

Degree=(1:numel(Z_local_mean))';

plot_sigma=plot(1e5*1000*(Z_local_mean+Z_local_std),'-r','LineWidth',1);
plot_sigma_shape=plot(1e5*1000*(Z_local_mean_shape+Z_local_std_shape),'--r','LineWidth',1);


Z_local_mean_std=Z_local_std/sqrt(NTapers);
Z_local_mean_std_shape=Z_local_std_shape/sqrt(NTapers);

plot(1e5*1000*(Z_local_mean-Z_local_mean_std),'-r','LineWidth',1);
plot(1e5*1000*(Z_local_mean_shape-Z_local_mean_std_shape),'--r','LineWidth',1);

% errorbar(Degree,Z_local_mean,Z_local_std,'LineWidth',3,'Color','r')


ylim(1e5*1000*[min(Z_global) max(Z_global)])

set(gca,'FontSize',20);

legend([plot_global, plot_global_shape,...
    plot_local_mean,plot_local_mean_shape,plot_sigma,plot_sigma_shape],...
{'global','global homo','local mean','local mean homo','1-sigma','1-sigma homo'});


drhodx=rho_mean./Z_local_mean_shape;
drhody=-rho_mean*Z_local_mean./(Z_local_mean_shape.^2);


density_local_std=sqrt((drhodx.^2).*(Z_local_mean_std.^2)+(drhody.^2).*(Z_local_mean_std_shape.^2))



%% Density plot

figure; hold on; box on;
set(gca,'FontSize',20);
density_global=Z_global./Z_global_shape*rho_mean;
density_local=Z_local_mean./Z_local_mean_shape*rho_mean;



plot(density_global,'-ob','LineWidth',3)


p=polyfit(3:18,density_global(3:18),1)

density_lin=polyval(p,3:18);

plot(3:18,density_lin,'--k','LineWidth',3)

plot(density_local,'-or','LineWidth',3)

% plot(density_local+density_local_std,'--k','LineWidth',2);
% plot(density_local-density_local_std,'--k','LineWidth',2);

legend({'Global','Localized in Rheasilvia','1-sigma'},'FontSize',20);

xlabel('Degree','FontSize',20);
ylabel('Effective density [kg/m^3]','FontSize',20);

xlim([2 20]);



density_local_mean=mean(density_local(6:14))



%% Localized correlation

Sgg_cross_shape_local_mean=mean(Sgg_cross_shape_local(:,Tapers),2);
Sgg_shape_local_mean=mean(Sgg_shape_local(:,Tapers),2);
Sgg_local_mean=mean(Sgg_local(:,Tapers),2);

R_local_mean=Sgg_cross_shape_local_mean./sqrt(Sgg_shape_local_mean.*Sgg_local_mean);

figure(Corr_fig)

plot_local_corr=plot(1:numel(R_local_mean),R_local_mean,'-or','LineWidth',5);


%% plot total consentration power

% V=1-V;

[loni,lati]=meshgrid(lor,lar);

% [F,lo,la,Plms]=plm2xyz([jk1 jk2 C{Tapers(1)}],1); F=F.^2*V(1);
[F,lo,la,Plms]=plm2xyz(lmcosi_window{Tapers(1)},1); F=F.^2*V2(1);

for i=2:numel(Tapers)
    %F=F+plm2xyz([jk1 jk2 C{Tapers(i)}],1,[],[],[],Plms).^2*V(Tapers(i));
    F=F+plm2xyz(lmcosi_window{Tapers(i)},1,[],[],[],Plms).^2*V2(Tapers(i));
    
end

AGUaxes

plotm(latc/180*pi,lonc/180*pi,'--k','LineWidth',4);

pcolorm(lati/180*pi,loni/180*pi,F);

cbar=colorbar;
ylabel(cbar,'Energy','FontSize',20)


%% Plotting Eigenfunctions

% Tapers=find((V2>0.6))

Vsort=-sort(-V2);

figure
set(gca,'FontSize',20);
plot(Vsort,'-o','MarkerSize',5);
xlabel('Taper number','FontSize',20);
ylabel('Concentration factor','FontSize',20);
xlim([1 36])

AGUaxes


h=pcolorm(lati/180*pi,loni/180*pi,r2(:,:,i));

plotm(latc/180*pi,lonc/180*pi,'--k','LineWidth',4);

title(['Concentration factor = ' num2str(V2(1))],'FontSize',20);

% plotm(Boundary(:,2)/180*pi,Boundary(:,1)/180*pi,'-k','LineWidth',3);

% lighting phong
% shading interp
vidObj = VideoWriter('Tapers.avi');
open(vidObj);


% Close the file.



for i=1:numel(Tapers)
    
    pause(1)
%     set(h,'CData',r(:,:,i).*ri);
    Tapers(i)=find(V2==Vsort(i),1);
    V2(Tapers(i))=-1;
    set(h,'CData',r2(:,:,Tapers(i)));    
    title(['Concentration factor = ' num2str(Vsort(i))],'FontSize',20);
    
currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);

end

close(vidObj);