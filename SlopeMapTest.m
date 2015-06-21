tic
%%

% ro_crust=2900;
% ro_mantle=3200;
% ro_core=6000;

load('DataForSlopes.mat');

ro_mantle_diff=ro_mantle-ro_crust;
ro_core_diff=ro_core-ro_mantle; 

Rref=265000;
omega=0.000326718;

%% Load Tri Shape Model
load SH-Tri_shape_model_uniform.mat

%% Load SH Shape Model

% MaxDegreeTopo=720;
% 
% load VestaHASTALAVESTAshape_sh1500.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% 
% %% Create a mesh
% Step=.2;
% 
% [lambdai,fii]=meshgrid(0:Step:360,90:-Step:-90);
% 
% lambdai=lambdai/180*pi;
% fii=fii/180*pi;
% 
% ri=plm2xyz(lmcosi_shape,Step);
% 
% [xi,yi,zi]=sph2cart(lambdai,fii,ri);
% 
% % MapRadialGrid(flipud(ri));
%  
% %  hold on;
% figure
% surf(xi,yi,zi,ri);
% StandardLight

%% Or load grided shape model

ShapeModelFileName='~/Dawn/Balmino/VestaTest/VestaPOSTVESTAshape_geoc_elev12m.grd';
% finfo = ncinfo(ShapeModelFileName)

[lambdai,fii,ri]=ReadGRD(ShapeModelFileName);

ri=ri*1000;
lambdai=lambdai/180*pi;
fii=fii/180*pi;

[xi,yi,zi]=sph2cart(lambdai,fii,ri);

% figure
% surf(xi,yi,zi,ri);
% StandardLight

%% Load SH Shape of core and mantle


% lmcosi_core_grav_new=data.lmcosi_core_grav_new;
% lmcosi_mantle_grav_new=data.lmcosi_mantle_grav_new
% mu_core_diff=data.mu_core_diff
% mu_mantle_diff=data.mu_mantle_diff;
% trisurf([Ver1' Ver2' Ver3'],x,y,z);
% alpha(0.5)
%  axis equal

%% Compute acceleration
% [U_p,ax_p,ay_p,az_p]=GravityPotentialTri2(x,y,z,Ver1,Ver2,Ver3,xi,yi,zi,ro_crust);
% 
% ax_p=-ax_p;
% ay_p=-ay_p;
% az_p=-az_p;
% 
% save('GravtityAccelerationTri02deg_2.mat','ax_p','ay_p','az_p','U_p');

%% Or load accelerations

data=load('GravtityAccelerationTri02deg_2.mat');

ax_p=data.ax_p;
ay_p=data.ay_p;
az_p=data.az_p;

clear data

%% Compute accelerations from mantle and core

[ax_core,ay_core,az_core]=GravityAcceleration(mu_core_diff,Rref,lmcosi_core_grav_new,xi,yi,zi);
[ax_mantle,ay_mantle,az_mantle]=GravityAcceleration(mu_mantle_diff,Rref,lmcosi_mantle_grav_new,xi,yi,zi);

ax_core=-ax_core;
ay_core=-ay_core;
az_core=-az_core;

ax_mantle=-ax_mantle;
ay_mantle=-ay_mantle;
az_mantle=-az_mantle;

a_mantle=sqrt(ax_mantle.^2+ay_mantle.^2+az_mantle.^2);
a_core=sqrt(ax_core.^2+ay_core.^2+az_core.^2);
a_p=sqrt(ax_p.^2+ay_p.^2+az_p.^2);

ac=omega*omega*sqrt(xi.^2+yi.^2);

ax=ax_p+ax_core+ax_mantle+omega*omega*sqrt(xi.^2+yi.^2).*cos(lambdai);
ay=ay_p+ay_core+ay_mantle+omega*omega*sqrt(xi.^2+yi.^2).*sin(lambdai);
az=az_p+az_core+az_mantle;

a=sqrt((ax.^2)+(ay.^2)+(az.^2));

axn=ax./a;
ayn=ay./a;
azn=az./a;

[nxi,nyi,nzi] = surfnorm(xi,yi,zi);

nxi=-nxi;
nyi=-nyi;
nzi=-nzi;

mean_nxi=(nxi(:,1)+nxi(:,end))/2;
mean_nyi=(nyi(:,1)+nyi(:,end))/2;
mean_nzi=(nzi(:,1)+nzi(:,end))/2;

mean_norm=sqrt(mean_nxi.^2+mean_nyi.^2+mean_nzi.^2);

mean_nxi=mean_nxi./mean_norm;
mean_nyi=mean_nyi./mean_norm;
mean_nzi=mean_nzi./mean_norm;

nxi(:,1)=mean_nxi;
nyi(:,1)=mean_nyi;
nzi(:,1)=mean_nzi;

nxi(:,end)=mean_nxi;
nyi(:,end)=mean_nyi;
nzi(:,end)=mean_nzi;

mean_nxi=mean(nxi(1,:));
mean_nyi=mean(nyi(1,:));
mean_nzi=mean(nzi(1,:));

nxi(1,:)=nxi(1,:)*0+mean_nxi;
nyi(1,:)=nyi(1,:)*0+mean_nyi;
nzi(1,:)=nzi(1,:)*0+mean_nzi;

mean_nxi=mean(nxi(end,:));
mean_nyi=mean(nyi(end,:));
mean_nzi=mean(nzi(end,:));

mean_norm=sqrt(mean_nxi.^2+mean_nyi.^2+mean_nzi.^2);

mean_nxi=mean_nxi/mean_norm;
mean_nyi=mean_nyi/mean_norm;
mean_nzi=mean_nzi/mean_norm;

nxi(end,:)=nxi(end,:)*0+mean_nxi;
nyi(end,:)=nyi(end,:)*0+mean_nyi;
nzi(end,:)=nzi(end,:)*0+mean_nzi;

% hold on;
% ToShow=fix(rand(1,100).*numel(xi));
% 
% 
% quiver3(xi(ToShow),yi(ToShow),zi(ToShow),nxi(ToShow),nyi(ToShow),nzi(ToShow));
% quiver3(xi(ToShow),yi(ToShow),zi(ToShow),ax_core(ToShow),ay_core(ToShow),az_core(ToShow),'r');
% quiver3(xi(ToShow),yi(ToShow),zi(ToShow),ax_mantle(ToShow),ay_mantle(ToShow),az_mantle(ToShow),'g');
% quiver3(xi(ToShow),yi(ToShow),zi(ToShow),ax_p(ToShow),ay_p(ToShow),az_p(ToShow),'b');
% quiver3(xi(ToShow),yi(ToShow),zi(ToShow),...
%     omega*omega*sqrt(xi(ToShow).^2+yi(ToShow).^2).*cos(lambdai(ToShow)),...
%     omega*omega*sqrt(xi(ToShow).^2+yi(ToShow).^2).*sin(lambdai(ToShow)),...
%     0*zi(ToShow),'k');
% 
% axis equal
% 
% hold on;
% quiver3(xi(:,1),yi(:,1),zi(:,1),nxi(:,1),nyi(:,1),nzi(:,1),'b');
% quiver3(xi(:,end),yi(:,end),zi(:,end),nxi(:,end),nyi(:,end),nzi(:,end),'r');
% axis equal;

ag=280000;
bg=226000;

[s_up,s_east,s_north]=GravityComponents(nxi,nyi,nzi,xi,yi,zi,ag,bg);

slope_az=atan2(s_north,s_east);



ang_diff=acos(nxi(:,1).*nxi(:,end) + nyi(:,1).*nyi(:,end) + nzi(:,1).*nzi(:,end));
slope=acos((axn.*nxi+ayn.*nyi+azn.*nzi))*180/pi;

% ds=slopeSPCmod-slopeSPG;
% hist(slope(:));

MapRadialGrid(slope);

MapRadialGrid(slope_az);


figure; hold on;
set(gca,'FontSize',20);
hist(log10(abs(ds(:))),50);
xlabel('log_{10}(\Delta slope) [deg]','FontSize',20);
ylabel('# of points','FontSize',20);
box on;

% MapRadialGrid(flipud(a_p)*1e5);
% MapRadialGrid(flipud(a_mantle)*1e5);
% MapRadialGrid(flipud(a_core)*1e5);
% MapRadialGrid(flipud(a)*1e5);

% slope in sh
% [lmcosi_slope,~]=xyz2plm(slope,2,'im',[],[],[]);
% slope_sh=plm2xyz(lmcosi_slope,1);
% MapRadialGrid(slope_sh);

%% curvature

% [gm, samc] = mcurvature_vec(xi,yi,zi);
% 
% MapRadialGrid(log(abs(gm)),-15,-5)
% 
% CurvGridFileName='~/Dawn/Balmino/VestaTest/VestaCurvGrid.xyz';
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,log(abs(gm)),CurvGridFileName);
% 
% % curv is sh
% 
% [lmcosi_curv]=xyz2plm(log(abs(gm)),8,'im',[],[],[]);
% 
% [curv_sh,lon,lat]=plm2xyz(lmcosi_curv,.2);
% 
% [lon,lat]=meshgrid(lon,lat);
% 
% MapRadialGrid(flipud(curv_sh));
% 
% CurvGridFileName_sh='~/Dawn/Balmino/VestaTest/VestaCurvGrid_sh.xyz';
% 
% WriteXYZ(lambdai*180/pi,fii*180/pi,curv_sh,CurvGridFileName_sh);

%% Write netCDF file

SlopeGridFileName='~/Dawn/Balmino/VestaTest/VestaSlopeGrid.xyz';

WriteXYZ(lambdai*180/pi,fii*180/pi,slope,SlopeGridFileName);


%% Print Histogram
% N_rand_points=1000000;
% [fi_rand,lambda_rand]=GenerateRandomSphericalCoord(N_rand_points);
% 
% 
% 
% % slope_col=slope(:);
% % lambda_col=lambdai(:);
% % fii_col=fii(:);
% 
% slope_rand=griddata(lambdai,fii,slope,lambda_rand,fi_rand,'nearest');
% 
% 
% slope_rand(slope_rand>45)=[];
% 
% 
% 
% slope_fig_hist = figure('XVisual',''); 
% axes1 = axes('Parent',slope_fig_hist,'FontSize',7);
% 
% 
% hold on;
% x_slope=0:1:40;
% 
% % hist(slope(slope<45),x_slope,50','b');
% n_elements = histc(slope_rand,x_slope);
% xlabel('Slope','FontSize',7);
% ylabel('[%]','FontSize',7);
% 
% xlim([0 40]);
% 
% ylabels = get(gca, 'YTickLabel');
% ylabels = linspace(0,100,length(ylabels));
% set(gca,'YTickLabel',ylabels);
% 
% set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
% set(gcf, 'PaperPositionMode','auto')
% 
% grid off;
% 
% box(axes1,'on');
% slope_elements = cumsum(n_elements);
% 
% % stairs(x_slope,slope_elements,'-r','LineWidth',2);
% hist(slope_rand,50);
%  
% NumberOfSlopes=numel(slope_rand); 
% y_max=max(max(get(gca,'YTick'))); 
% max_ytick=y_max/NumberOfSlopes*100;
% ylabels = get(gca, 'YTickLabel');
% ylabels = linspace(0,100,length(ylabels));
% set(gca,'YTickLabel',(ylabels*max_ytick/100));
% set(gca,'YTickLabel',fix(10*ylabels*max_ytick/100)/10);
% 
% 
% print(slope_fig_hist, '-dpsc', 'SlopeHist_c.eps');

%% plot 3D

% h=surf(xi,yi,zi,slope);
% 
% 
% lighting phong
% shading interp
% axis equal tight off
% set(gcf,'Color','k');
% box off
% hold on;
% light_handle=light('Style','infinite');
% set(h,'BackFaceLighting','lit');
% % colormap gray;
% material([0 0.9 0]);

 toc