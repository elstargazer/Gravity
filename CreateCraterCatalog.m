ccc

%% Input parameters

a=280900;
c=226200;
Ncr=20;
R=265000;

MaxDeg=20;
Resolution=0.05;

latlim=[-90  -50   0  50  90];
latcen=[  -90 -25  25  90];
lonlim=[0  60  120  180  240  300  360];
loncen=[ 30  90  150  210  270  330];

%% Load Shape model
ShapeModelFileName1='~/Dawn/Balmino/VestaTest/Vesta20140513shape_geoc_elev_3m.grd';
ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR64_3m.grd';
% ShapeModelFileName1='~/Dawn/Balmino/VestaTest/VestaDTMCut.grd';

[lambdai1,fii1,ri1]=ReadGRD(ShapeModelFileName1); ri1=ri1*1000;
[xi1,yi1,zi1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);

[lambdai2,fii2,ri2]=ReadGRD(ShapeModelFileName2); 
[xi2,yi2,zi2]=sph2cart(lambdai2/180*pi,fii2/180*pi,ri2);

lmcosi=xyz2plm(ri1,MaxDeg,'im');
ri1_sh=plm2xyz(lmcosi,Resolution);
% e=Eccentricity(a,c);
% [~,~,H1]=XYZ2BLH(xi1,yi1,zi1,a,e);
H1=ri1-ri1_sh;
H2=ri2-ri1_sh;

%% Show map

latin=3;
lonin=1;

Hm=(H1-min(H1(:)))/(max(H1(:))-min(H1(:)));

fig=figure('Color','w','Position',[1 1 1400 1400]);

axm=axesm('stereo','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],'FontSize',24,'Grid','on','MLabelParallel',...
    'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
    , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,'FontColor',[0 0 0]...
    ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0],...
    'FLatLimit',[latlim(latin) latlim(latin+1)],'MapLonLimit',[lonlim(lonin) lonlim(lonin+1)],...
    'MapLatLimit',[latlim(latin) latlim(latin+1)],'FLonLimit',[lonlim(lonin) lonlim(lonin+1)],...
    'origin',[latcen(latin) loncen(lonin)]);

%% Pick craters

s=3;
surfm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi...
    ,H1(1:s:end,1:s:end),Hm(1:s:end,1:s:end)); shading interp;
light('Position',[-1 1 1],'Style','infinite');

latcc=zeros(1,Ncr);
latco=zeros(1,Ncr);
loncc=zeros(1,Ncr);
lonco=zeros(1,Ncr);
rc=zeros(1,Ncr);

for i=1:Ncr
    
    [latc, lonc] = inputm(2);
    
    latcc(i)=latc(1)*180/pi;
    loncc(i)=lonc(1)*180/pi;    
    latco(i)=latc(2)*180/pi;
    lonco(i)=lonc(2)*180/pi;
       
    rc(i)=distance(latc(1)*180/pi,lonc(1)*180/pi,...
        latc(2)*180/pi,lonc(2)*180/pi);
    
    [latsc,lonsc] = scircle2(latc(1)*180/pi,lonc(1)*180/pi,latc(2)*180/pi,lonc(2)*180/pi);
    
    plot3m(latsc/180*pi,lonsc/180*pi,2,'-k','LineWidth',3);  
    drawnow    
end

% AGUaxes;

%% Run crater depth estimation

[depth1,rc_e1]=CraterDepth(fii1,lambdai1,H1,latcc,loncc,rc);
[depth2,rc_e2]=CraterDepth(fii2,lambdai2,H2,latcc,loncc,rc);

%% Crater Depth VS Crater diameter

% figure; hold on;
% set(gca,'FontSize',20);
% 
% plot(2*rc,depth1/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
% plot(2*rc,depth2/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% 
% % plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
% 
% xlabel('Crater diameter [km]','FontSize',20);
% ylabel('Crater depth [km]','FontSize',20);
% legend({'SPC','SPG'},'FontSize',20);


figure; hold on;
set(gca,'FontSize',20);

plot(2*rc,(depth1-depth2)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');

% plot(2*rc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

xlabel('Crater diameter [km]','FontSize',20);
ylabel('Crater depth difference[km]','FontSize',20);
% legend({'SPC','SPG'},'FontSize',20);


%% Crater Depth VS latitude

figure; hold on;
set(gca,'FontSize',20);

plot(latcc,depth1/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(latcc,depth2/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

% plot(latcc,(depth2-depth1)/1000,'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

xlabel('Latitude [deg]','FontSize',20);
ylabel('Crater depth [km]','FontSize',20);
legend({'SPC','SPG'},'FontSize',20);








