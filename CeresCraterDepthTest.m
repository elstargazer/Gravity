ccc

%% Input parameters

aref=482000;
cref=446000;
R=470000;

MaxDeg=20;
Resolution=1/16;

%% Load Shape model
shape_folder1='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename1='SHAPE_SPC150828_512.bds';

shape_folder2='/Users/antonermakov/Dawn/CeresShapeModel/SPG/Survey/';
shape_filename2='global.bds';

[~,shapename2,~] = fileparts(shape_filename2) ;
full_filename2 = [shape_folder2 shape_filename2];
[~,shapename1,~] = fileparts(shape_filename1) ;
full_filename1 = [shape_folder1 shape_filename1];

% crater data base
crater_base_filename='~/Dawn/ceres_totalcratercount_equi_Survey.dat';

%% Read crater calalog

[latcc,loncc,diam]=ReadCraterCatalog(crater_base_filename);
Ncr=numel(diam);
log10diam=log10(diam);
radcc=(diam/2)/(R/1000)*180/pi;

%%
AGUaxes;
for i=1:numel(radcc)    
    [latsc,lonsc] = scircle1(latcc(i),loncc(i),radcc(i));    
    plot3m(latsc,lonsc,2,'-k','LineWidth',2);   
end

[x,y,z]=ReadSPC(full_filename1,0.5,'grid');
[loni,lati,ri] = cart2sph(x,y,z);

[B,L,H] = XYZ2BLH(x,y,z,aref/1000,Eccentricity(aref,cref));

s=1;
surfm(lati(1:s:end,1:s:end)*180/pi,loni(1:s:end,1:s:end)*180/pi...
    ,ri(1:s:end,1:s:end),H(1:s:end,1:s:end)); shading interp;
light('Position',[-1 1 1],'Style','infinite');

cbar=colorbar('FontSize',20);
ylabel(cbar,'Elevation [m]','FontSize',20);

%% Run crater depth estimation

MaxDiam=50;
MinDiam=10;

mask=(diam<MaxDiam) & (diam>MinDiam);

latcc_t=latcc(mask);
loncc_t=loncc(mask);
radcc_t=radcc(mask);
diam_t=diam(mask);

Ncr=numel(radcc_t);

[depth1,depth1_std,rc_e1]=CraterDepthDSK(full_filename1,latcc_t,loncc_t,radcc_t);
[depth2,depth2_std,rc_e2]=CraterDepthDSK(full_filename2,latcc_t,loncc_t,radcc_t);

depth2 = depth2/1000;
depth2_std = depth2_std/1000;
rc_e2 = rc_e2 / 1000;

depth1_relstd=depth1_std./depth1;
depth2_relstd=depth2_std./depth2;

dd=(depth1-depth2);

dd_std=sqrt(depth1_std.^2+depth2_std.^2);

%% PLOTTING THINGS
%% Crater Depth VS Crater diameter

figure; hold on;
set(gca,'FontSize',20);

plot(diam(mask),depth1,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(diam(mask),depth2,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

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

[binCenters,binMean,binStandardDev]=rebin(diam(mask),(depth1-depth2),1);

plot(binCenters,binMean)
h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',2);
% errorbar_tick(h,1000);
grid on;
title('SPC - SPG','FontSize',20);

set(gca,'XScale','log');

%% Crater Depth VS latitude 

figure; hold on;
set(gca,'FontSize',20);

plot(latcc_t,depth1,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(latcc_t,depth2,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

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

[binCenters,binMean,binStandardDev]=rebin(latcc_t,(depth1-depth2),2.5);

plot(binCenters,binMean)
h=errorbar(binCenters,binMean,binStandardDev,'LineWidth',2);
% errorbar_tick(h,1000);
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

%% Plot Crater depth difference on map
% H1=H2;
% 
% s=4;
% 
% minH1=min(H1(1:s:end,1:s:end)); minH1=min(minH1);
% maxH1=max(H1(1:s:end,1:s:end)); maxH1=max(maxH1);
% 
% Hm=(H1(1:s:end,1:s:end)-minH1)/(maxH1-minH1);
% 
% mindd=min(dd_nn)/2;
% maxdd=max(dd_nn)/2;
% 
% % fig=figure('Color','w','Position',[1 1 1400 1400]);
% 
% fig=figure('Color','k','Position',[1 1 1400 1400]); hold on;
% 
% title('SPC - SPG','FontSize',20);
% 
% sb2=subplot('Position',[0.05 0.04 0.9 0.05]);
% axis off
% 
% colormap jet
% cbar=colorbar('FontSize',20,'Location','South');
% xlabel(cbar,'Crater depth difference [km]','FontSize',20,'Color','w');
% caxis([mindd maxdd])
% 
% 
% sb1=subplot('Position',[0.05 0.12 0.9 0.6]);
% 
% ax=axesm('mollweid','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],'FontSize',24,'Grid','on','MLabelParallel',...
%     'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
%     , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,'FontColor',[1 1 1]...
%     ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0],'FEdgeColor',[1 1 1],'Frame','on');
% 
% axis off;

% for i=1:Ncr    
%     [latsc,lonsc] = scircle1(latcc(i),loncc(i),radcc(i));    
%     plot3m(latsc/180*pi,lonsc/180*pi,2,'-k','LineWidth',2);   
%     i/Ncr
% end

% scatterm(latcc_nn/180*pi,loncc_nn/180*pi,30,dd_nn,'o','filled');
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Crater depth difference [km]','FontSize',20);
% caxis([-1 1]/5);
% 
% cc2=jet(128);
% 
% colin=fix((dd_nn-mindd)/(maxdd-mindd)*127)+1;
% 
% colin(colin<1)=1;
% colin(colin>128)=128;
% 
% 
% for i=1:numel(latcc_nn)    
%     [latsc,lonsc] = scircle1(latcc_nn(i),loncc_nn(i),radcc_nn(i));   
%     if ~isnan(colin(i))
%         plot3m(latsc/180*pi,lonsc/180*pi,2,'-','LineWidth',2,'color',cc2(colin(i),:));  
% %         plot3m(latsc/180*pi,lonsc/180*pi,2,'-','LineWidth',2,'color','k');
%     end
%     i/Ncr
% end
% 
% freezeColors
% cbfreeze(cbar)

% % cc2(end,:) = [0;0;0]; %black 
% colormap(cc2); %activate it.
% caxis([mindd maxdd])

% colormap grey
% 
% sf=surfm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi...
%     ,H1(1:s:end,1:s:end),Hm); shading interp;
% light('Position',[-1 1 1],'Style','infinite');

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



