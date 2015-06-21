ccc

%% Input Paramerers;
GravityFileName='/Users/antonermakov/GRAIL/Gravity_Models/jggrx_0660b_sha.tab.txt';
% GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

MaxDegreeTopo=1800;
MaxDegreeGrav=500;
% Rref_gravity=293000;
Rref_gravity=1738000;
rhomean=3344;
Resolution=.1;
MaxTopoPower=4;
L=20;
MinConcentration=0.93;
NTess=2;
circle_rad=10; 
MaxDegreeExp=7;
rho_mean=3457.5;
GoodDegrees=2+L:MaxDegreeGrav-L;

%% Load gravity model

[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

% load lmcosi_calc_new.mat
% lmcosi_grav=lmcosi_calc_new;
% lmcosi_grav=TruncateGravityModel(lmcosi_grav,MaxDegreeGrav,0);

% lmcosi_grav(3,3)=0;

lmcosi_fa=plm2pot(lmcosi_grav,Rref_gravity,mu,Rref,3,'nothing');

% plotting geoid
% lmcosi_geoid=plm2pot(lmcosi_grav,Rref_gravity,mu,Rref,0,'nothing');
% 
% [r_geoid,lon_geoid,lat_geoid]=plm2xyz(lmcosi_geoid,0.1);
% [lon_geoid,lat_geoid]=meshgrid(lon_geoid,lat_geoid);
% 
% [x_geoid,y_geoid,z_geoid]=sph2cart(lon_geoid/180*pi,lat_geoid/180*pi,r_geoid*1000+1737000);
% 
% figure('color','k'); hold on;
% h=surf(x_geoid,y_geoid,z_geoid,r_geoid); shading interp;
% 
% lighting phong
% shading interp
% axis equal tight off
% box off
% hold on;
% light_handle=light('Style','infinite');
% set(h,'BackFaceLighting','lit');
% colormap jet;
% material([0 0.9 0]);
% view(180,0)

%% Load topography model

% load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=load('/Users/antonermakov/GRAIL/Topography/MoonTopo2600p.shape.sh');

lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% lmcosi_grav(4,3)=0;

[ri,~,~,~]=plm2xyz(lmcosi_shape,Resolution);
MeanRadius=lmcosi_shape(1,3);
ri=ri-MeanRadius;
MapRadialGrid(flipud(ri));

%% Homogeneous gravity



lmcosi_grav_shape=TopoSH2GravitySH(ri+MeanRadius,mu,rho_mean,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);
[U,~,~,~]=plm2xyz(lmcosi_grav_shape,Resolution,[],[],[]);

lmcosi_grav_shape(1,:)=[];
lmcosi_fa_shape=plm2pot(lmcosi_grav_shape,Rref_gravity,mu,Rref,3,'nothing');

%% Computing free-air anomaly

% [fa,~,~,~]=plm2xyz(lmcosi_fa,Resolution,[],[],[]);
% [fa_shape,~,~,~]=plm2xyz(lmcosi_fa_shape,Resolution,[],[],[]);

%% Plotting Free-air anomaly

% [fa_fig,fa_ax]=MapRadialGrid(flipud(fa)*1e5);
% MapRadialGrid(flipud(ri))
 
%% Global spectral power and admittance

% lmcosi_shape(:,3:4)=lmcosi_shape(:,3:4)./lmcosi_shape(1,3);
% lmcosi_fa(:,3:4)=lmcosi_fa(:,3:4)./lmcosi_fa(1,3);

[fa,~,~,~]=plm2xyz(lmcosi_fa,Resolution,[],[],[]);
[ri,~,~,~]=plm2xyz(lmcosi_shape,Resolution);

Z_global=SphericalHarmonicAdmittance(lmcosi_fa,lmcosi_shape);
% Z_global_shape=SphericalHarmonicAdmittance(lmcosi_fa_shape,lmcosi_shape);

Sgt_global=CrossPower(lmcosi_fa,lmcosi_shape);
Stt_global=CrossPower(lmcosi_shape,lmcosi_shape);
Sgg_global=CrossPower(lmcosi_fa,lmcosi_fa);

C_global=Sgt_global./sqrt(Stt_global.*Sgg_global(1:size(Stt_global,2)));

% Sgg_global_shape=CrossPower(lmcosi_fa_shape,lmcosi_fa_shape);
% Sgg_global_cross_shape=CrossPower(lmcosi_fa,lmcosi_fa_shape);

figure('Position',[1 1 1000 700]/1.3); hold on;
set(gca,'FontSize',20);

degree=1:MaxDegreeGrav;
plot(degree,Z_global*1e5*1000,'-or','LineWidth',2,'MarkerSize',3);

xlabel('Degree []','FontSize',20);
ylabel('Admittance [mGal/km]','FontSize',20);


haxes1=gca;
haxes1_pos = get(haxes1,'Position'); % store position of first axes

haxes2 = axes('Position',haxes1_pos,...
              'YAxisLocation','right',...
              'Color','none','FontSize',20,'XTick',[]);
         
% plot(haxes2,fii1(:,1),rlonav,'-k','LineWidth',3)
line(degree,C_global,'Parent',haxes2,'Color','b','LineWidth',4)
ylabel(haxes2,'Correlation []','FontSize',20);

ylim([-1 1])
box on;

% plot(degree,Z_global_shape,'-or','LineWidth',4,'MarkerSize',4);

g=1.6;
nu=0.25;
rhocrust=2600;
rhomantle=3220;
E=40e9;
d=20000;
tl=2*d;
n=1:100;

% Z_flex=AiryAdmittance(Rref,1.6,0.25,rhomean,rhomantle,rhocrust,E,d,n,tl);
% plot(n,Z_flex,'-or','LineWidth',3,'MarkerSize',3);

x0=[4e10 20000 2600];
lb=[1e9 5000 2300];
ub=[1e11 100000 3600];

[x_fit,resnorm] = lsqcurvefit(@AiryAdmittanceFun,x0,degree,Z_global*1e8,lb,ub);
Z_fit=AiryAdmittanceFun(x_fit,degree);
plot(degree,Z_fit,'-or','LineWidth',1,'MarkerSize',3);

%% Global correlation
% R_global=Sgg_global_cross_shape./sqrt(Sgg_global_shape.*Sgg_global);

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

% fii=fii(1:10);
Npoints=numel(fii);
% MaxDegreeExp=fix(0.5*(-3+sqrt(1+8*Npoints)))

%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,0,0);

figure; hold on;
set(gca,'FontSize',20);

plot(-sort(-V2),'-ok','LineWidth',3);

xlabel('Taper number','FontSize',20);
ylabel('Concentration factor [] ','FontSize',20);
box on;
grid on;
xlim([1 numel(V2)]);

s=size(G2);

for i=1:s(2)       
    lmcosi_window_basic{i}=glm2lmcosi(G2,i);    
end

Tapers=find(V2>MinConcentration);
NumberOfTapers=numel(Tapers);
lmcosi_window_basic=lmcosi_window_basic(Tapers);
NTapers=sum(V2>MinConcentration);

Z_t=zeros(numel(GoodDegrees),numel(fii));
Z_t_std=zeros(numel(GoodDegrees),numel(fii));
C_t=zeros(numel(GoodDegrees),numel(fii));
C_t_std=zeros(numel(GoodDegrees),numel(fii));
Z_tm=zeros(numel(GoodDegrees),NTapers*numel(fii));
C_tm=zeros(numel(GoodDegrees),NTapers*numel(fii));

progressbar(0);

for j=1:numel(fii)

% patch coordinates
    lat_center=fii(j);
    lon_center=lambdai(j);

    colat_center=90-lat_center;
    alp=0;
    bta=lat_center-90; 
    gam=-lon_center;
    
%     [fig_cm,ax_cm]=AGUaxes;
%     set(fig_cm,'Position',[1 1 1000*1.6 500*1.6]);
%     vidObj = VideoWriter('Tapers_Moon.avi');
%     vidObj.FrameRate=5;
%     vidObj.Quality=100;
%     open(vidObj);

    for i=1:NumberOfTapers
        
        [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,bta,gam,'dlmb');    
        [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},Resolution);    
    
        [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,MaxDegreeGrav,'im');    
        [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,MaxDegreeGrav,'im');
        
        
        %plotting tapers
        [lori,lari]=meshgrid(lor/180*pi,lar/180*pi);
        
        
%         AGUaxes
%         [latcirc,loncirc]=scircle1(lat_center,lon_center,circle_rad);
%         pcolorm(lari,lori,r(:,:,i)); shading interp;
%         plotm(latcirc/180*pi,loncirc/180*pi,'-k','LineWidth',3)
%         caxis([-1.5 1.5]);
%         title(['Taper # ' num2str(i) ' '],'FontSize',24);
        
%         lmcosi_shape_w{i}(:,3:4)=lmcosi_shape_w{i}(:,3:4)./lmcosi_shape_w{i}(1,3);
%         lmcosi_fa_w{i}(:,3:4)=lmcosi_fa_w{i}(:,3:4)./lmcosi_fa_w{i}(1,3);
         
        Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
        Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    
        Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
            
%          currFrame = getframe(gcf);
%          writeVideo(vidObj,currFrame);

    end   
    
%     close(vidObj);

%% Localized admittance plot

    N0=(L+1)*(circle_rad/180);
    Sgt_local_mean=mean(Sgt_local,2);
    Stt_local_mean=mean(Stt_local,2);
    Sgg_local_mean=mean(Sgg_local,2);
    
    Z_local_mean=Sgt_local_mean./Stt_local_mean;
    Z_local_std=std(Sgt_local./Stt_local,0,2);
    Z_local_mean_std=Z_local_std/sqrt(NTapers);
    Degree=(1:numel(Z_local_mean))';
    Z_t(:,j)=Z_local_mean(GoodDegrees);
    Z_t_std(:,j)=Z_local_mean_std(GoodDegrees);
    
    Z_tm(:,(j-1)*NTapers+1:j*NTapers)=Sgt_local(GoodDegrees,:)./Stt_local(GoodDegrees,:);
    C_tm(:,(j-1)*NTapers+1:j*NTapers)=Sgt_local(GoodDegrees,:)./...
        sqrt(Stt_local(GoodDegrees,:).*Sgg_local(GoodDegrees,:));   
    
    C_local_mean=Sgt_local_mean./sqrt(Stt_local_mean.*Sgg_local_mean);
    C_local_std=std(Sgt_local./sqrt(Stt_local.*Sgg_local),0,2);
    C_local_mean_std=C_local_std/sqrt(NTapers);
    C_t(:,j)=C_local_mean(GoodDegrees);
    C_t_std(:,j)=C_local_mean_std(GoodDegrees);
 
    progressbar(j/numel(fii));
    j/numel(fii)

end
progressbar(1);

%% Fitting

x0=[4e10 20000 2600];
lb=[1e9 5000 2300];
ub=[1e11 100000 3600];

progressbar(0);
clear x_fit resnorm

for j=1:numel(fii)

    [x_fit(j,:),resnorm(j)] = lsqcurvefit(@AiryAdmittanceFun,x0,GoodDegrees',1e8*Z_t(:,j),lb,ub);
    Z_t_fit(:,j)=AiryAdmittanceFun(x_fit(j,:),GoodDegrees');
    progressbar(j/numel(fii));

end

progressbar(1);

%% Plotting

AGUaxes

for i=1:numel(fii)    
    [latc,lonc] = scircle1(fii(i),lambdai(i),circle_rad);
    
%     if ((fii(i)<10) & (fii(i)>-10) & (lambdai(i) <10) & (lambdai(i)>-10))
%      plotm(latc/180*pi,lonc/180*pi,'-r','LineWidth',4); 
      plotm(fii(i)/180*pi,lambdai(i)/180*pi,'ob','LineWidth',1,'MarkerFaceColor','b');
%     end
end
% plotm(fii/180*pi,lambdai/180*pi,'ob','MarkerFaceColor','b','MarkerSize',10);

%% Making a movie

figs=figure('Position',[1 1 1700 1300]);
hold on

s1=subplot('Position',[0.08 0.05 0.84 0.5]);
set(gca,'FontSize',20);
s2=subplot('Position',[0.08 0.65 0.84 0.30]);
set(gca,'FontSize',20);

subplot(s1)

ax=axesm('mollweid','frame','on','FontSize',24,'Grid','on','MLabelParallel',...
    'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel',...
     'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',2,'FontColor',[0 0 0],...
     'FontSize',24,'GAltitude',Inf,'Geoid',[1 0]);

s=size(ri);
fit=linspace(-90,90,s(1))/180*pi;
lambdat=linspace(0,360,s(2))/180*pi;
[lambdati,fiti]=meshgrid(lambdat,fit);

surfm(fiti,lambdati,flipud(ri));
shading interp

[latc,lonc] = scircle1(fii(1),lambdai(1),circle_rad);
plot_sc=plotm(latc/180*pi,lonc/180*pi,'-r','LineWidth',4);  
subplot(s2);

hold on;

plot_Z=plot(GoodDegrees,1e8*Z_t(:,1),'-ok','MarkerSize',4,'LineWidth',3);
%plot_Z_fit=plot(GoodDegrees,Z_t_fit(:,1),'-or','MarkerSize',3,'LineWidth',1);
h=errorbar(GoodDegrees,1e8*Z_t(:,1),1e8*Z_t_std(:,1),'LineWidth',1,'Color','k');
errorbar_tick(h,1000);

xlim([0 200]);
ylim([-200 200])

xlabel('Degree','FontSize',20);
ylabel('Admittance [mGal/km]','FontSize',20);

vidObj = VideoWriter('MoonLocalAdm3.avi');
vidObj.FrameRate=5;
vidObj.Quality=100;

open(vidObj);


[~,index]=sortrows(floor([fii,lambdai]*100)/100,[1 2])

for i=1:numel(fii)  
    
    unplot(1);
    
         
    [latc,lonc] = scircle1(fii(index(i)),lambdai(index(i)),circle_rad);
    
    subplot(s1)    
    delete(plot_sc)    
    plot_sc=plotm(latc/180*pi,lonc/180*pi,'-r','LineWidth',4);
    plotm(fii(index(i))/180*pi,lambdai(index(i))/180*pi,'o','MarkerSize',10,'MarkerFaceColor','k');      
    set(plot_Z,'YData',1e8*Z_t(:,i));
    subplot(s2)
    h=errorbar(GoodDegrees,1e8*Z_t(:,i),1e8*Z_t_std(:,i),'LineWidth',1,'Color','k');
    errorbar_tick(h,1000);
    
    %unplot(1);
    %set(plot_Z_fit,'YData',Z_t_fit(:,i));
    
%     drawnow;
%     pause(0.1)
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);    
    i
      
end

close(vidObj);

%% Plot all Z
figure; hold on;
set(gca,'FontSize',20);

plot(GoodDegrees,Z_t*1e8,'LineWidth',3)

xlabel('SH Degree','FontSize',20);
ylabel('Admittance [mGal/km] ','FontSize',20);
box on;
grid on;


%% Plot all Corr

figure; hold on;
set(gca,'FontSize',20);

plot(GoodDegrees,C_t,'LineWidth',3)

xlabel('SH Degree','FontSize',20);
ylabel('Correlation [] ','FontSize',20);
box on;
grid on;

%% Clusters

NClucters=4;

% Admittance clusters

[idx,ctrs,sumd,D] = kmeans(Z_t',NClucters,'distance','sqEuclidean','onlinephase','on');
cc=jet(NClucters);
Dm=log10(min(D'));
[madm,ixs]=sort(mean(ctrs,2));
[~,ixs2]=sort(ixs);
cc2=cc(ixs2,:);

for i=1:NClucters
    ctrs_std(:,i)=std(Z_t(:,idx==i),0,2);
end

% Correlation clusters

[idx_c,ctrs_c,sumd_c,D_c] = kmeans(C_t',NClucters,'distance','sqEuclidean','onlinephase','on');
cc_c=jet(NClucters);
Dm_c=log10(min(D_c'));
[madm_c,ixs_c]=sort(mean(ctrs_c,2));
[~,ixs2_c]=sort(ixs_c);
cc2_c=cc(ixs2_c,:);

for i=1:NClucters
    ctrs_c_std(:,i)=std(C_t(:,idx_c==i),0,2);
end


cc3=jet(128);
idc=fix(((Dm-min(Dm))/(max(Dm)-min(Dm))*127))+1;

% hist(Dm(idx==3),30)
for i=1:NClucters
    Dm_med(i)=median(Dm(idx==i));
end

%% Plot centers admittance
plot_z=figure; hold on;
set(gca,'FontSize',20);


for i=1:NClucters    
    plot(GoodDegrees,1e8*ctrs(i,:),'-o','MarkerSize',3,'LineWidth',2,'Color',cc2(i,:));
%     h=errorbar(GoodDegrees,1e8*ctrs(i,:),1e8*ctrs_std(:,i),'LineWidth',1,'Color',cc2(i,:));
%     errorbar_tick(h,1000);
end

xlabel('SH Degree','FontSize',20);
ylabel('Admittance [mGal/km] ','FontSize',20);
box on;
grid on;

legend({'1','2','3','4'},'FontSize',20);

%% Plot centers correlation
plot_c=figure; hold on;
set(gca,'FontSize',20);

for i=1:NClucters    
    plot(GoodDegrees,ctrs_c(i,:),'-o','MarkerSize',3,'LineWidth',2,'Color',cc2_c(i,:));
    h=errorbar(GoodDegrees,ctrs_c(i,:),ctrs_c_std(:,i),'LineWidth',1,'Color',cc2_c(i,:));
    errorbar_tick(h,1000);
end

xlabel('SH Degree','FontSize',20);
ylabel('Correlation []','FontSize',20);
box on;
grid on;


% ctrs_c_std

% x0=[4e10 20000 30000];
% lb=[1e9 0 0];
% ub=[1e11 100000 100000];

x0=[4e10 20000 2600];
lb=[1e9 0 2000];
ub=[1e11 100000 4000];

% (R,g,nu,rhomean,rhomantle,rhocrust,E,d,n,tl)

Z_f=AiryAdmittanceFun(x0,GoodDegrees);

% clear x_f_fit,Z_f_fit;

for j=1:NClucters

    [x_f_fit(j,:),~] = lsqcurvefit(@AiryAdmittanceFun,x0,GoodDegrees',1e8*ctrs(j,:)',lb,ub);
    Z_f_fit(:,j)=AiryAdmittanceFun(x_f_fit(j,:),GoodDegrees');

end


% plot(GoodDegrees,Z_f_fit,'k--','LineWidth',3);

figure(plot_z);

for i=1:NClucters    
    plot(GoodDegrees,Z_f_fit(:,i)','-.k','MarkerSize',1,'LineWidth',3,'MarkerFaceColor',cc2(i,:));
end

text_legend=cell(NClucters,1);
for i=1:NClucters    
    text_legend{i}=num2str(ixs2(i));
end

legend(text_legend,'FontSize',20);

x_f_fit(:,1)=x_f_fit(:,1)/1e9;


in=fopen('Z_fit_param.txt','w');


for i=1:NClucters
    fprintf(in,'Ncl = %d, E = %E [GPa], d = %f [m], rho = %f [km/m^3]\n',[i x_f_fit(i,:)]);
end

fclose(in);

%% Plot Voronoi 
AGUaxes;
scatterm(fii/180*pi,lambdai/180*pi,100,idx,'filled');

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

alpha(0.5);

axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);

%% Mapping hexagons

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    plotm(fi_b,lambda_b,'-r','LineWidth',2);
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end

%% Mapping minimum distance

AGUaxes; hold on
% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');

in_Dm=fix((Dm-min(Dm))./(max(Dm)-min(Dm))*127)+1;

cc_Dm=jet(128)

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    fill3m(fi_b,lambda_b,1,cc_Dm(in_Dm(k),:),'EdgeColor','none');
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end

cbar=colorbar('FontSize',20);
ylabel(cbar,'Deviation from cluster [log_{10} (s^{-2})]','FontSize',20)
caxis([min(Dm) max(Dm)])

%% Mapping clusters admittance
AGUaxes; hold on
% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');

ncl = size(clmap,1);

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    fillm(fi_b,lambda_b,0,cc2(idx(k),:),'EdgeColor','none');
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end

%% Mapping clusters correlation
AGUaxes; hold on
% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');

ncl = size(clmap,1);

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    fillm(fi_b,lambda_b,0,cc2_c(idx_c(k),:),'EdgeColor','none');
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end


%% Record Z classes

lambda=(0:0.2:360);
fi=(-90:0.2:90);
[lambdap,fip]=meshgrid(lambda,fi);
z_class=zeros(size(fip));

progressbar(0);

parfor i=1:numel(fip)    
    [dist,az] = distance(fip(i),lambdap(i),fii,lambdai);    
    j=(dist==min(dist));    
    z_class(i)=idx(find(j,1));    
%     progressbar(i/numel(fip));
end

progressbar(1);

z_class2=z_class;

for i=1:NClucters    
    z_class2(z_class==ixs(i))=i;    
end

AGUaxes
pcolorm(fip/180*pi,lambdap/180*pi,z_class2); shading interp;
WriteXYZ(lambdap,fip,z_class2,'/Users/antonermakov/GRAIL/Topography/Z_class2.xyz');

%% Record C classes

lambda=(0:.2:360);
fi=(-90:.2:90);
[lambdap,fip]=meshgrid(lambda,fi);
c_class=zeros(size(fip));

progressbar(0);

parfor i=1:numel(fip)    
    [dist,az] = distance(fip(i),lambdap(i),fii,lambdai);    
    j=(dist==min(dist));    
    c_class(i)=idx_c(find(j,1));    
    progressbar(i/numel(fip));
end

progressbar(1);

c_class2=c_class;

for i=1:NClucters    
    c_class2(c_class==ixs_c(i))=i;    
end

AGUaxes
pcolorm(fip/180*pi,lambdap/180*pi,c_class2); shading interp;
WriteXYZ(lambdap,fip,c_class2,'/Users/antonermakov/GRAIL/Topography/C_class2.xyz');

%% Movie admittance on surface
Lmax=20;

vidObj = VideoWriter('MoonAdm_comp.avi','Motion JPEG AVI');
vidObj.FrameRate=5;
vidObj.Quality=30;
open(vidObj);

lmcosi_adm=xyz2plm(Z_t(1,:),Lmax,'irr',fii',lambdai');    
[Zi,lon,lat]=plm2xyz(lmcosi_adm,1);
    
[lon,lat]=meshgrid(lon,lat);

lon=lon/180*pi;
lat=lat/180*pi;

[fig_cm,ax_cm]=AGUaxes;
set(fig_cm,'Position',[1 1 1000*1.6 500*1.6]);
% PlotContours
i=1;
plot_zi=pcolorm(lat,lon,Zi*1e8); shading interp; 
cbar=colorbar('FontSize',20);
ylabel(cbar,'Admittance [mGal/km] ','FontSize',20);
title(['Degree = ' num2str(GoodDegrees(i)) ' '],'FontSize',20);
caxis([-200 200]);

currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);

for i=2:numel(GoodDegrees)
    
    lmcosi_adm=xyz2plm(Z_t(i,:),Lmax,'irr',fii',lambdai');    
    Zi=plm2xyz(lmcosi_adm,1);
    
    set(plot_zi,'CData',Zi*1e8);
    title(['Degree = ' num2str(GoodDegrees(i)) ' '],'FontSize',20);     

    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end

close(vidObj);

%% Movie correlation on surface
Lmax=20;

vidObj = VideoWriter('MoonCorr_comp.avi','Motion JPEG AVI');
vidObj.FrameRate=5;
vidObj.Quality=30;
open(vidObj);

lmcosi_corr=xyz2plm(C_t(1,:),Lmax,'irr',fii',lambdai');    
[Ci,lon,lat]=plm2xyz(lmcosi_corr,1);
    
[lon,lat]=meshgrid(lon,lat);

lon=lon/180*pi;
lat=lat/180*pi;

[fig_cm,ax_cm]=AGUaxes;
set(fig_cm,'Position',[1 1 1000*1.6 500*1.6]);
i=1;
plot_ci=pcolorm(lat,lon,Ci); shading interp; 
cbar=colorbar('FontSize',20);
ylabel(cbar,'Correlation [] ','FontSize',20);
title(['Degree = ' num2str(GoodDegrees(i)) ' '],'FontSize',20);
caxis([-1 1]);

currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);

for i=2:numel(GoodDegrees)
    
    lmcosi_adm=xyz2plm(C_t(i,:),Lmax,'irr',fii',lambdai');    
    Ci=plm2xyz(lmcosi_adm,1);
    
    set(plot_ci,'CData',Ci);
    title(['Degree = ' num2str(GoodDegrees(i)) ' '],'FontSize',20);
     

    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end     

close(vidObj);


