ccc

%% Input parameters

a=280900;
c=226200;
R=265000;

MaxDeg=20;
Resolution=1/64;

%% Load Shape model
ShapeModelFileName1='/Users/ermakov/Dawn/ShapeModel/Gaskell/20140513/Vesta20140513shape_geoc_elev.grd';
ShapeModelFileName2='/Users/ermakov/Dawn/ShapeModel/DLR/DLR_SHAPE_2014_06_19/VestaShapeDLR64.grd';

% ShapeModelFileName1='~/Dawn/Balmino/VestaTest/Vesta20130522shape_geoc_elev_3m_gridline.grd';
% ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR_3m.grd';
 
   
[lambdai1,fii1,ri1]=ReadGRD(ShapeModelFileName1); ri1=ri1*1000;
[xi1,yi1,zi1]=sph2cart(lambdai1(1:16:end,1:16:end)/180*pi,...
    fii1(1:16:end,1:16:end)/180*pi,ri1(1:16:end,1:16:end));


[~,~,ri2]=ReadGRD(ShapeModelFileName2); ri2=flipud(ri2);
% [xi2,yi2,zi2]=sph2cart(lambdai2/180*pi,fii2/180*pi,ri2);

lmcosi=xyz2plm(ri1(1:16:end,1:16:end),MaxDeg,'im');
ri1_sh=plm2xyz(lmcosi,Resolution);
% e=Eccentricity(a,c);
% [~,~,H1]=XYZ2BLH(xi1,yi1,zi1,a,e);
H1=ri1-ri1_sh;
H2=ri2-ri1_sh;

clear ri1_sh

figure
pcolor(H1(1:16:end,1:16:end)); shading interp
figure
pcolor(H2(1:16:end,1:16:end)); shading interp

%% plotting shape models

% fig=figure('Color','w','Position',[1 1 1400 1400]); hold on;
% 
% ax=axesm('mollweid','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],'FontSize',24,'Grid','on','MLabelParallel',...
%     'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
%     , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,'FontColor',[0 0 0]...
%     ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0],'FEdgeColor',[0 0 0],'Frame','on');
% 
% axis off;
% 
% s=16;
% 
% H1=ri1;
% 
% minH1=min(H1(1:s:end,1:s:end)); minH1=min(minH1);
% maxH1=max(H1(1:s:end,1:s:end)); maxH1=max(maxH1);
% 
% Hm=(H1(1:s:end,1:s:end)-minH1)/(maxH1-minH1);
% 
% sf=surfm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi...
%     ,H1(1:s:end,1:s:end),Hm); shading interp;
% light('Position',[-1 1 1],'Style','infinite');
% 
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation [m]','FontSize',20);

%% Read crater calalog

filename='~/Dawn/marchi_GCCv2_bis.dat';

[latcc,loncc,diam]=ReadCraterCatalog(filename);
Ncr=numel(diam);
log10diam=log10(diam);
radcc=(diam/2)/(R/1000)*180/pi;

%% Plot craters
% % H1=H2;
% Hm=(H1-min(H1(:)))/(max(H1(:))-min(H1(:)));
% 
% % fig=figure('Color','w','Position',[1 1 1400 1400]);
% 
% AGUaxes;
% 
% for i=1:Ncr    
%     [latsc,lonsc] = scircle1(latcc(i),loncc(i),radcc(i));    
%     plot3m(latsc/180*pi,lonsc/180*pi,2,'-k','LineWidth',2);   
%     i/Ncr
% end
% 
% s=1;
% surfm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi...
%     ,H1(1:s:end,1:s:end),Hm(1:s:end,1:s:end)); shading interp;
% light('Position',[-1 1 1],'Style','infinite');
% 
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation [m]','FontSize',20);

%% Run crater depth estimation

MaxDiam=100;
MinDiam=0;

mask=(diam<MaxDiam) & (diam>MinDiam);

latcc_t=latcc(mask);
loncc_t=loncc(mask);
radcc_t=radcc(mask);
diam_t=diam(mask);

Ncr=numel(radcc_t);

[depth1,depth1_std,rc_e1]=CraterDepth(fii1,lambdai1,H1,latcc_t,loncc_t,radcc_t);
[depth2,depth2_std,rc_e2]=CraterDepth(fii1,lambdai1,H2,latcc_t,loncc_t,radcc_t);

depth1_relstd=depth1_std./depth1;
depth2_relstd=depth2_std./depth2;

dd=(depth1-depth2)/1000;

dd_std=sqrt(depth1_std.^2+depth2_std.^2);

%% PLOTTING THINGS
%% Crater Depth VS Crater diameter

figure; hold on;
set(gca,'FontSize',20);

plot(diam(mask),depth1/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(diam(mask),depth2/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

% h=errorbar(diam(mask),depth1/1000,depth1_std/1000,'LineWidth',2);
% errorbar_tick(h,1000);
% 
% h=errorbar(diam(mask),depth2/1000,depth2_std/1000,'LineWidth',2);
% errorbar_tick(h,1000);

% plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

xlabel('Crater diameter [km]','FontSize',20);
ylabel('Crater depth [km]','FontSize',20);
legend({'SPC','SPG'},'FontSize',20);

% plot difference
figure; hold on;
set(gca,'FontSize',20);
plot(diam(mask),dd,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
xlabel('Crater diameter [km]','FontSize',20);
ylabel('Crater depth difference [km]','FontSize',20);

[binCenters,binMean,binStandardDev]=rebin(diam(mask),(depth1-depth2)/1000,1);

plot(binCenters,binMean)
h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',2);
errorbar_tick(h,1000);
grid on;
title('SPC - SPG','FontSize',20);

set(gca,'XScale','log');

%% Crater Depth VS latitude 

figure; hold on;
set(gca,'FontSize',20);

plot(latcc_t,depth1/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(latcc_t,depth2/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

% plot(latcc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

xlabel('Latitude [deg]','FontSize',20);
ylabel('Crater depth [km]','FontSize',20);
legend({'SPC','SPG'},'FontSize',20);
xlim([-90 90]);

% plot difference

figure; hold on;
set(gca,'FontSize',20);
plot(latcc_t,dd,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
xlabel('Latitude [deg]','FontSize',20);
ylabel('Crater depth difference [km]','FontSize',20);
xlim([-90 90]);
set(gca,'XTick',-90:30:90);

[binCenters,binMean,binStandardDev]=rebin(latcc_t,(depth1-depth2)/1000,2.5);

plot(binCenters,binMean)
h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',2);
errorbar_tick(h,1000);
grid on;
title('SPC - SPG','FontSize',20);

%% Crater Depth VS latitude VS diameter

figure; hold on;
set(gca,'FontSize',20);
scatter((diam(mask)),latcc_t,40,dd,'filled');
xlabel('Diameter [km]','FontSize',20);
ylabel('Latitude [deg]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Depth difference [km]')
ylim([-90 90]);
set(gca,'YTick',-90:30:90);
set(gca,'XScale','log')
caxis([-0.75 0.75]);
xlim([MinDiam MaxDiam]);
title('SPC - SPG','FontSize',20);

%% 2D rebin

xdata=log10(diam(mask));
ydata=latcc_t;
zdata=dd;
zdata_std=dd_std;

sx=0.05;
sy=2.5;

[XbinCenters,YbinCenters,binMean,binStandardDev]...
    =rebin2D( xdata,ydata,zdata,zdata_std,sx,sy );

[YbinCenters,XbinCenters]=meshgrid(YbinCenters,XbinCenters);

% mean plot
figure; hold on;
pcolor(XbinCenters,YbinCenters,binMean); shading flat;

xlabel('Diameter [log_{10} (km)]','FontSize',20);
ylabel('Latitude [deg]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Mean depth difference [km]')

ylim([-90 90]);
set(gca,'YTick',-90:30:90);

xlim(log10([MinDiam MaxDiam]));
title('SPC - SPG','FontSize',20);
caxis([-2 2]);
box on;

% std plot
figure; hold on;
pcolor(XbinCenters,YbinCenters,binStandardDev); shading flat;

xlabel('Diameter [log_{10} (km)]','FontSize',20);
ylabel('Latitude [deg]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Depth difference standard deviation [km]')

ylim([-90 90]);
set(gca,'YTick',-90:30:90);

xlim(log10([MinDiam MaxDiam]));
title('SPC - SPG','FontSize',20);
caxis([0 1]);
box on


%% Histogram

figure('Color','k'); hold on;
title('SPC - SPG','FontSize',20,'Color','w');
set(gca,'FontSize',20,'Color','w');


[n1, xout1] = hist(dd,50);
bar(xout1,n1,'b'); 

[n2, xout2] = hist(dd(latcc_t<40),50);
bar(xout2,n2,'r');

alpha(0.5);

xlabel('Crater depth difference [km]','FontSize',20,'Color','w');

box on;

legend({'All craters','Below 50N'},'FontSize',20,'Color','w');
ylabel('Number of craters','FontSize',20,'Color','w');

xlim([-0.75 0.75]);

%% Surface gravity at craters
% 
% MaxDeg=100;
% omega=1617.3331237/180*pi/86400;
% 
% % load tri shape model
% load 'VestaTriShape4.mat'
% ri=sqrt(xi.*xi+yi.*yi+zi.*zi);
% 
% 
% % shape model in SH
% lmcosi=xyz2plm(ri1(1:16:end,1:16:end),MaxDeg,'im');
% 
% % compute radius vector at craters
% 
% latcc_mod=-latcc_t;
% loncc_mod=loncc_t+180;
% % loncc_mod(loncc_mod<0)=loncc_mod(loncc_mod<0)+180;
% 
% rcr=plm2xyz(lmcosi,latcc_mod,loncc_mod); %rcr=flipud(rcr);
% [xcr,ycr,zcr]=sph2cart(loncc_t/180*pi,latcc_t/180*pi,rcr'+1000);
% 
% tic
% [ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV,xcr,ycr,zcr,3457*ones(size(FV,1),1)');
% ax=ax-omega*omega*sqrt(xcr.^2+ycr.^2).*cos(loncc_t/180*pi);
% ay=ay-omega*omega*sqrt(xcr.^2+ycr.^2).*sin(loncc_t/180*pi);
% atot=sqrt(ax.*ax+ay.*ay+az.*az);
% toc
% 
% %% Crater Depth VS Crater diameter vs g
% 
% MaxRelSTD=0.7;
% MaxDiam2=30;
% 
% figure; hold on;
% set(gca,'FontSize',20);
% 
% % scatter(diam(mask),depth2/1000,30,atot,'o','filled');
% % h=errorbar(diam(mask),depth1/1000,depth1_std/1000,'LineWidth',2);
% % errorbar_tick(h,1000);% 
% % h=errorbar(diam(mask),depth2/1000,depth2_std/1000,'LineWidth',2);
% % errorbar_tick(h,1000);
% % plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% 
% mask2=(~isnan(depth1)) & (~isnan(depth1_std)) & ...
%     (depth1_relstd<MaxRelSTD) & (diam_t<MaxDiam2) & (latcc_t<90);
% 
% depth1_nn=depth1(mask2);
% depth1_relstd_nn=depth1_relstd(mask2);
% diam_nn=diam_t(mask2);
% atot_nn=atot(mask2);
% latcc_nn=latcc_t(mask2);
% loncc_nn=loncc_t(mask2);
% radcc_nn=radcc_t(mask2);
% dd_nn=dd(mask2);
% 
% p=polyfit(diam_nn,depth1_nn,1);
% 
% diam_lin=0:30;
% 
% depth1_fit=polyval(p,diam_nn);
% depth1_lin=polyval(p,diam_lin);
% 
% scatter(diam_nn,depth1_nn/1000,30,atot_nn,'o','filled');
% 
% % plot(diam_nn,depth1_fit/1000,'.k');
% plot(diam_lin,depth1_lin/1000,'--k','LineWidth',4);
% 
% 
% xlabel('Crater diameter [km]','FontSize',20);
% ylabel('Crater depth [km]','FontSize',20);
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Gravitational Acceleration [m/sec^2] ','FontSize',20);
% ylim([0 7]);
% xlim([0 MaxDiam2]);
% box on;
% 
% % y=LinFun(a,x)
% 
% f=@(x,xdata) x*xdata;
% slope = lsqcurvefit(f,p(1),diam_nn,depth1_nn);
% depth1_lin0=f(slope,diam_lin);
% depth1_fit0=f(slope,diam_nn);
% 
% plot(diam_lin,depth1_lin0/1000,'-.k','LineWidth',4);
% 
% %% fits with acceleration bins
% 
% ai=0.21:0.01:0.26;
% cc=jet(numel(ai)-1);
% 
% figure; hold on;
% set(gca,'FontSize',20);
% 
% for i=1:numel(ai)-1
%     
%     mask3=((atot_nn>ai(i)) & (atot_nn<ai(i+1)));  
%     
% %     pi=polyfit(diam_nn(mask3),depth1_nn(mask3),1);    
% %     depth1_lin=polyval(pi,diam_lin);
% %     plot(diam_lin,depth1_lin/1000,'--','LineWidth',4,'Color',cc(i,:));
% 
%     pause(0.5);
%     drawnow;
%     
%     slopei(i) = lsqcurvefit(f,p(1),diam_nn(mask3),depth1_nn(mask3));   
%     depth1_lin0=f(slopei(i),diam_lin);     
%     plot(diam_lin,depth1_lin0/1000,'-.','LineWidth',4,'Color',cc(i,:));    
%     
% end
% 
% xlabel('Crater diameter [km]','FontSize',20);
% ylabel('Crater depth [km]','FontSize',20);
% 
% ylim([0 7]);
% xlim([0 MaxDiam2]);
% box on;
% 
% figure; hold on;
% set(gca,'FontSize',20);
% plot(slopei)
% box on;
% 
% %% plot map with acceleration
% 
% figure; hold on;
% set(gca,'FontSize',20);
% scatter(loncc_nn,latcc_nn,30,atot_nn,'o','filled');
% 
% xlabel('Longitude [deg]','FontSize',20);
% ylabel('Latitude [deg]','FontSize',20);
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Gravitational Acceleration [m/sec^2] ','FontSize',20);
% xlim([-180 180]);
% ylim([-90 90]);
% box on;
% 
% % plot map with deviation
% figure; hold on;
% set(gca,'FontSize',20);
% scatter(loncc_nn,latcc_nn,30,(depth1_nn-depth1_fit0)/1000,'o','filled');
% 
% xlabel('Longitude [deg]','FontSize',20);
% ylabel('Latitude [deg]','FontSize',20);
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Deviation from linear law [km] ','FontSize',20);
% xlim([-180 180]);
% ylim([-90 90]);
% caxis([-1 1]);
% box on;
% 
% %% Deviation plot
% figure; hold on
% set(gca,'FontSize',20);
% 
% plot(atot_nn,(depth1_nn-depth1_fit)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% xlabel('Gravitational Acceleration [m/sec^{2}]','FontSize',20);
% ylabel('Deviation from linear law [km]','FontSize',20);
% 
% grid on;
% 
% [binCenters,binMean,binStandardDev]=rebin(atot_nn,(depth1_nn-depth1_fit)/1000,0.01);
% 
% plot(binCenters,binMean)
% h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',4);
% errorbar_tick(h,1000);
% box on;
% 
% %% depth to diameters ratio
% 
% MaxRelSTD=0.7;
% MaxDiam2=30;
% 
% figure; hold on;
% set(gca,'FontSize',20);
% 
% % scatter(diam(mask),depth2/1000,30,atot,'o','filled');
% % h=errorbar(diam(mask),depth1/1000,depth1_std/1000,'LineWidth',2);
% % errorbar_tick(h,1000);% 
% % h=errorbar(diam(mask),depth2/1000,depth2_std/1000,'LineWidth',2);
% % errorbar_tick(h,1000);
% % plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% 
% mask2=(~isnan(depth1)) & (~isnan(depth1_std)) & ...
%     (depth1_relstd<MaxRelSTD) & (diam_t<MaxDiam2) & (latcc_t<90);
% 
% depth1_nn=depth1(mask2);
% depth1_relstd_nn=depth1_relstd(mask2);
% diam_nn=diam_t(mask2);
% atot_nn=atot(mask2);
% latcc_nn=latcc_t(mask2);
% loncc_nn=loncc_t(mask2);
% 
% dtd_nn=depth1_nn./diam_nn;
% 
% p=polyfit(diam_nn,dtd_nn,1);
% 
% diam_lin=0:30;
% 
% dtd_nn_fit=polyval(p,diam_nn);
% dtd_nn_lin=polyval(p,diam_lin);
% 
% inrand=fix(rand(1,numel(diam_nn))*numel(diam_nn))+1;
% 
% 
% scatter(diam_nn(inrand),dtd_nn(inrand)/1000,30,atot_nn(inrand),'o','filled');
% % scatter(diam_nn,dtd_nn/1000,30,atot_nn,'o','filled');
% 
% % plot(diam_nn,depth1_fit/1000,'.k');
% plot(diam_lin,dtd_nn_lin/1000,'--k','LineWidth',4);
% 
% 
% xlabel('Crater depth to diameter ratio []','FontSize',20);
% ylabel('Crater depth [km]','FontSize',20);
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Gravitational Acceleration [m/sec^2] ','FontSize',20);
% ylim([0 0.5]);
% xlim([0 MaxDiam2]);
% 
% caxis([0.21 0.245])
% box on;
% 
% %% fits with acceleration bins depth to diameter
% 
% ai=min(atot_nn):0.01:max(atot_nn);
% cc=jet(numel(ai)-1);
% 
% % figure; hold on;
% % set(gca,'FontSize',20);
% 
% for i=1:numel(ai)-1
%     
%     mask3=((atot_nn>ai(i)) & (atot_nn<ai(i+1)));  
%     
%     p_i=polyfit(diam_nn(mask3),dtd_nn(mask3),1);    
%     dtd_lini=polyval(p_i,diam_lin);
%     plot(diam_lin,dtd_lini/1000,'--','LineWidth',4,'Color',cc(i,:));
% 
% end
% 
% ylabel('Crater depth to diameter ratio []','FontSize',20);
% xlabel('Crater depth [km]','FontSize',20);
% 
% %% DTD deviation plot
% 
% figure; hold on
% set(gca,'FontSize',20);
% 
% plot(atot_nn,(dtd_nn-dtd_nn_fit)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% xlabel('Gravitational Acceleration [m/sec^{2}]','FontSize',20);
% ylabel('DTD deviation from linear law []','FontSize',20);
% grid on;
% 
% [binCenters,binMean,binStandardDev]=rebin(atot_nn,(dtd_nn-dtd_nn_fit)/1000,0.005);
% plot(binCenters,binMean)
% h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',4);
% errorbar_tick(h,1000);
% box on;
% 
% xlim([0.21 0.26]);
% ylim([-0.3 0.3]);
% 
% %% 
% 
% figure; hold on
% set(gca,'FontSize',20);
% 
% plot(atot_nn,(dtd_nn)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% p=polyfit(atot_nn,(dtd_nn),1);
% 
% alin=0.20:0.001:0.26;
% dtd_nn_lina=polyval(p,alin);
% 
% plot(alin,(dtd_nn_lina)/1000,'--b','LineWidth',4);
% xlabel('Gravitational Acceleration [m/sec^{2}]','FontSize',20);
% ylabel('DTD []','FontSize',20);
% grid on;
% 
%% Plot Crater depth difference on map
% H1=H2;

s=4;

minH1=min(H1(1:s:end,1:s:end)); minH1=min(minH1);
maxH1=max(H1(1:s:end,1:s:end)); maxH1=max(maxH1);

Hm=(H1(1:s:end,1:s:end)-minH1)/(maxH1-minH1);

mindd=min(dd_nn)/2;
maxdd=max(dd_nn)/2;

% fig=figure('Color','w','Position',[1 1 1400 1400]);

fig=figure('Color','k','Position',[1 1 1400 1400]); hold on;

title('SPC - SPG','FontSize',20);

sb2=subplot('Position',[0.05 0.04 0.9 0.05]);
axis off

colormap jet
cbar=colorbar('FontSize',20,'Location','South');
xlabel(cbar,'Crater depth difference [km]','FontSize',20,'Color','w');
caxis([mindd maxdd])


sb1=subplot('Position',[0.05 0.12 0.9 0.6]);

ax=axesm('mollweid','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],'FontSize',24,'Grid','on','MLabelParallel',...
    'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
    , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,'FontColor',[1 1 1]...
    ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0],'FEdgeColor',[1 1 1],'Frame','on');

axis off;

% for i=1:Ncr    
%     [latsc,lonsc] = scircle1(latcc(i),loncc(i),radcc(i));    
%     plot3m(latsc/180*pi,lonsc/180*pi,2,'-k','LineWidth',2);   
%     i/Ncr
% end

% scatterm(latcc_nn/180*pi,loncc_nn/180*pi,30,dd_nn,'o','filled');
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Crater depth difference [km]','FontSize',20);
% caxis([-1 1]/5);

cc2=jet(128);

colin=fix((dd_nn-mindd)/(maxdd-mindd)*127)+1;

colin(colin<1)=1;
colin(colin>128)=128;


for i=1:numel(latcc_nn)    
    [latsc,lonsc] = scircle1(latcc_nn(i),loncc_nn(i),radcc_nn(i));   
    if ~isnan(colin(i))
        plot3m(latsc/180*pi,lonsc/180*pi,2,'-','LineWidth',2,'color',cc2(colin(i),:));  
%         plot3m(latsc/180*pi,lonsc/180*pi,2,'-','LineWidth',2,'color','k');
    end
    i/Ncr
end

freezeColors
cbfreeze(cbar)

% % cc2(end,:) = [0;0;0]; %black 
% colormap(cc2); %activate it.
% caxis([mindd maxdd])

colormap grey

sf=surfm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi...
    ,H1(1:s:end,1:s:end),Hm); shading interp;
light('Position',[-1 1 1],'Style','infinite');

% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation [m]','FontSize',20);

% cdata=zeros(size(Hm,1),size(Hm,2),3)+0.65;
% set(sf,'CData',cdata);

% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation [m]','FontSize',20);
% 
% colormap jet
% 
% latbin=-90:10:90;
% 
% for i=1:numel(latbin)-1       
%     stdi(i)=std(dd_nn(~isnan(dd_nn) & (latcc_nn<latbin(i+1)) & (latcc_nn>latbin(i)))) ;     
% end
% 
% %% STD of difference vs latitude
% 
% figure; hold on
% set(gca,'FontSize',20);
% 
% plot(latbin(1:end-1)+5,stdi,'-ok','LineWidth',4);
% 
% xlabel('Latitude [deg]','FontSize',20);
% ylabel('STD of difference [km]','FontSize',20);
% xlim([-90,90]);
% set(gca,'XTick',-90:30:90);
% box on;



