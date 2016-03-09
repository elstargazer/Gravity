ccc

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';

%% input parameters
% SPC shape
shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150828_512.bds';
filename_grav = '/Users/antonermakov/Dawn/CeresGravityModel/CERES18B01/JGC18B01.sha';

% SPG shape
% shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPG/Survey/';
% shape_filename='global.bds';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

% GM    = 62.6253e9;
G     = 6.67384e-11;
% Rref  = 476000;
aref  = 481000;
cref  = 446000;
step  = 0.1;
r1    = 470000;
T     = 9.073859324514187; % DLR
Npts  = 100;

MaxDegreeTopo = 100;
MaxDegreeGrav = 12;
MaxTopoPower  = 4;

%% get shape model in SH

[x_grid,y_grid,z_grid]=ReadSPC(full_filename,step,'grid');

x_grid=x_grid*1000;
y_grid=y_grid*1000;
z_grid=z_grid*1000;

V = Mesh2Volume(x_grid,y_grid,z_grid);

[lon_grid,lat_grid,r_grid]=cart2sph(x_grid,y_grid,z_grid);
% reference ellipsoid surface
[xref,yref,zref]=TriEllRadVec(lat_grid,lon_grid,aref,aref,cref,'xyz');

eccref=Eccentricity(aref,cref);
[B,L,H]=XYZ2BLH(x_grid,y_grid,z_grid,aref,eccref);

lmcosi_t = xyz2plm(flipud(r_grid'),MaxDegreeTopo);

%% get gravity model in SH

[lmcosi_g,Rref,GM,GM_std]=ReadGRAILGravityModel(filename_grav);
lmcosi_g = [0 0 1 0 0 0; lmcosi_g];

lmcosi_g = TruncateGravityModel(lmcosi_g,MaxDegreeGrav,1);
lmcosi_g(4,3) = lmcosi_g(4,3);
% GM = GM*1e9;
M=GM/G;
rhomean = M/V;

J2obs = -lmcosi_g(4,3);
lmcosi_gt1_ell = SHRotationalEllipsoid(481000,446000,2,Rref); 

% Plot gravity spectrum with error spectrum
h_grav_spec = figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

set(gca,'YScale','log');
set(gca,'XTick',2:MaxDegreeGrav);

[sdl_grav,l_sdl] = plm2spec(lmcosi_g);
[sdl_grav_std,l_sdl_std] = plm2spec([lmcosi_g(:,1:2) lmcosi_g(:,5:6)]);

plot(l_sdl,sdl_grav,'or-','LineWidth',3,'MarkerSize',5);
plot(l_sdl_std,sdl_grav_std,'ob-','LineWidth',3,'MarkerSize',5);

xlim([2 MaxDegreeGrav])

ylabel('Gravity PSD','FontSize',fntsize);
xlabel('Degree','FontSize',fntsize);

%% plot topography

% AGUaxes;
% pcolorm(lat_grid*180/pi,lon_grid*180/pi,H);

%% Plot topography with respect to equipotential surface

H_eq = Height2Equipotential(full_filename,lat_grid,lon_grid,GM,Rref,lmcosi_g,T);
% 
% AGUaxes;
% pcolorm(lat_grid*180/pi,lon_grid*180/pi,H_eq/1000);
% cbar = colorbar('FontSize',fntsize);
% ylabel(cbar,'Height above equipotential [km]','FontSize',20);
% caxis([-6 6]);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,H_eq/1000,'H_eq.dat');

%% hydrostatic gravity

% Grid of core radii and densities
% rcore=2500:500:470000;
% rhocoreg=2000:100:4000;

r2=linspace(10000,470000,Npts);
r2 = r2(1:end-1);
r2_add = linspace(r2(end),r1,10);
r2 = [r2 r2_add(2:end)];
rho2=linspace(rhomean,6000,Npts);

[rho2i,r2i]=meshgrid(rho2,r2);

rho1i=-(3*M-4*pi*(r2i.^3).*rho2i)./(4*pi*(r2i.^3)-4*pi*(r1^3));
rho1i(rho1i<0)=NaN;

M2 = 4/3*pi.*(r2i.^3).*(rho2i-rho1i);
M1 = 4/3*pi.*(r1.^3).*(rho1i);

% Compute hydrostatic flattening factors
[f2i,f1i]=HydrostaticStateExact2lGrid(r1,r2i,T,rho1i,rho2i);

J2hi=RadFlat2J2(r1,r2i,f1i,f2i,rho1i,rho2i,Rref);

fig_todel = figure;
CJhyd = contour(r2i,rho2i,J2hi,[J2obs J2obs]);
% CJhyd  = contourc(r2i(:,1),rho2i(1,:),J2hi,[J2obs J2obs]);
% plot(CJhyd(1,2:end),CJhyd(2,2:end),'-or');
% plot(CJhyda(1,2:end),CJhyda(2,2:end),'-ob');
close(fig_todel);

r2_Jh   = CJhyd(1,2:end-1);
rho2_Jh = CJhyd(2,2:end-1);

rho1_Jh  = griddata(r2i,rho2i,rho1i,r2_Jh,rho2_Jh,'cubic');
M2_Jh    = griddata(r2i,rho2i,M2,r2_Jh,rho2_Jh,'cubic');
fp2_Jh   = griddata(r2i,rho2i,f2i,r2_Jh,rho2_Jh,'cubic');
fp1_Jh   = griddata(r2i,rho2i,f1i,r2_Jh,rho2_Jh,'cubic');

save('2layer_solution.mat','rho1_Jh','rho2_Jh','r1','r2_Jh');

M1_Jh = M - M2_Jh;

fig_shell=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

plot(rho1_Jh,(r1-r2_Jh)/1000,'-b','LineWidth',3);

rho1_lin = 900:10:1900;

cond_nan = ~isnan(rho1_Jh);
st_lin = interp1(rho1_Jh(cond_nan),(r1-r2_Jh(cond_nan))/1000,rho1_lin,'cubic');
rho2_lin = interp1(rho1_Jh(cond_nan),rho2_Jh(cond_nan),rho1_lin,'cubic');
r2_lin = interp1(rho1_Jh(cond_nan),r2_Jh(cond_nan),rho1_lin,'cubic');

ans = interp1(rho1_Jh(cond_nan),rho2_Jh(cond_nan), 1361.9    ,'cubic');

% figure; hold on;
% plot(rho1_Jh(cond_nan),rho2_Jh(cond_nan),'-or');
% plot(rho1_lin,rho2_lin,'-ob');

xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);

% gi = ginput(1);

gi(1) = 1400;

ind = find(abs(rho1_Jh - gi(1)) == min(abs(rho1_Jh - gi(1))));
plot(rho1_Jh(ind),(r1-r2_Jh(ind))/1000,'or','MarkerSize',10);

PrintWhite(fig_shell,[fig_folder 'Fig_shell_pick.jpg']);

% we have r1, r2_Jh, rho1_Jh, rho2_Jh = family of solutions for J2

% open Ryan2layer.fig
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on;grid on; box on;
% 
% plot((r2_Jh(cond_nan))/1000,rho1_Jh(cond_nan),'-om')
% plot((r2_Jh(cond_nan))/1000,rho2_Jh(cond_nan),'-oc')
% 
% plot(r2_lin/1000,rho1_lin,'-ko')
% plot(r2_lin/1000,rho2_lin,'-ko')
% 
% xlabel('Core radius  [km]','FontSize',fntsize);
% ylabel('Density [kg/m^3]','FontSize',fntsize);
% 
% legend({'Shell','Core'},'FontSize',fntsize_sm);
% 
% in_2lmodel = fopen('2LayerModelJ2Grid.txt','w');
% fprintf(in_2lmodel,'rho1 (kg/m^3), rho2 (kg/m^3), r2 (km)\n');
% fprintf(in_2lmodel,'%6.2f, %6.2f, %6.2f\n',[rho1_lin; rho2_lin; r2_lin/1000]);
% fclose(in_2lmodel);

%% plot gravity

% compute just gravity 
[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_g,xref,yref,zref);
[g_up,g_east,g_north]=GravityComponents(ax,ay,az,xref,yref,zref,aref,cref);

% computing free-air anomaly
lmcosi_fa = lmcosi_g;
% lmcosi_fa(4,3) = 0;
% lmcosi_fa(1,3) = 0;

[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_fa,xref,yref,zref);
[g_up_fa,g_east_fa,g_north_fa]=GravityComponents(ax,ay,az,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,g_up_fa*1e5,'FA.dat');

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,g_up_fa*1e5); shading interp;
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Free-air anomaly [mGal]','FontSize',20);

PrintWhite([fig_folder 'Fig_FA.jpg']);

%% compute gravity from shape

% lmcosi_gt1=TopoSH2GravitySH(flipud(r_grid'),GM,rhomean,Rref,...
%     MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_gt1=Topo2Grav(flipud(r_grid'),Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

% lmcosi_gt1 = lmcosi_gtr;
% lmcosi_gt1(7:end,:) = [];
% lmcosi_gt1(1,3) = 1;

%% model gravity
% compute gravity from shape
[ax_gt1,ay_gt1,az_gt1]=GravityAcceleration(...
    GM,Rref,lmcosi_gt1,xref,yref,zref);
[g_up_gt1,g_east_gt1,g_north_gt1]=GravityComponents(...
    ax_gt1,ay_gt1,az_gt1,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,g_up_gt1*1e5,'GT.dat');

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,g_up_gt1*1e5)

figure;
surf(x_grid,y_grid,z_grid,g_up_gt1);
StandardLight;
axis equal;

a_Jh = zeros(size(r2_Jh));
c_Jh = a_Jh;

%% compute Bouguer anomaly

% compute gravity from core
[a_Jh,~,c_Jh]=fr2abc(r2_Jh(ind),fp2_Jh(ind),0);
lmcosi_gt2=SHRotationalEllipsoid(a_Jh,c_Jh,MaxDegreeGrav,Rref);
w = [M1_Jh(ind)/M M2_Jh(ind)/M];
    
lmcosi_gt = WeightSumExpansion(w,{lmcosi_gt1,lmcosi_gt2});

MaxDegreeBouguer=min([lmcosi_g(end,1) lmcosi_gt(end,1)]);
lmcosi_ba = lmcosi_g;
lmcosi_ba(:,3:4) = lmcosi_ba(:,3:4) - lmcosi_gt(:,3:4);

[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_ba,xref,yref,zref);
[gba_up,gba_east,gba_north]=GravityComponents(...
    ax,ay,az,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gba_up*1e5,'BA.dat');

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,gba_up*1e5); shading interp;
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Bouguer anomaly [mGal]','FontSize',20);

BA = gba_up*1e5;

image_BA = (BA - min(BA(:)))/...
    (max(BA(:)) - min(BA(:)))*(2^16-1);

imwrite(uint16(image_BA'),['Ceres_BA_BW_' ...
    num2str(min(BA(:))) '_' num2str(max(BA(:))) '.png']);

PrintWhite([fig_folder 'Fig_BA.jpg']);

%% Compute isostatic anomaly

D_comp = r1 - r2_Jh(ind);
 
% lmcosi_gtisos_lin = Topo2IsosGrav(...
%     flipud(r_grid'),Rref,D_comp,rho1_Jh(ind),rho2_Jh(ind),...
%     rhomean,MaxDegreeTopo,MaxDegreeGrav,1);
% 
% lmcosi_gtisos = Topo2IsosGrav(...
%     flipud(r_grid'),Rref,D_comp,rho1_Jh(ind),rho2_Jh(ind),...
%     rhomean,MaxDegreeTopo,MaxDegreeGrav,2);
% 
% MaxDegreeIsos=min([lmcosi_g(end,1) lmcosi_gtisos(end,1)]);
% lmcosi_isos        = lmcosi_g;
% lmcosi_isos(:,3:4) = lmcosi_isos(:,3:4) - lmcosi_gtisos(:,3:4);
% 
% [ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_isos,xref,yref,zref);
% [gisos_up,~,~]=GravityComponents(...
%     ax,ay,az,xref,yref,zref,aref,cref);
% 
% WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gisos_up*1e5,'ISOS.dat');
% 
% AGUaxes;
% pcolorm(lat_grid*180/pi,lon_grid*180/pi,gisos_up*1e5); shading interp;
% colorbar('FontSize',fntsize);
% ylabel(cbar,'Isostatic anomaly [mGal]','FontSize',20);
% 
% PrintWhite([fig_folder 'Fig_ISOS.jpg']);

% Isos_coef = ((r1./r2_Jh).^2).*rho1_Jh./(rho2_Jh-rho1_Jh);
% t = H_eq*Isos_coef(ind);

t = FindCrustalRoot(r1,D_comp,H_eq,rho1_Jh(ind),rho2_Jh(ind)-rho1_Jh(ind));

r2_grid = TriEllRadVec(lat_grid,lon_grid,a_Jh,a_Jh,c_Jh,'rad');
r2_grid = r2_grid - t;

% r2_grid = r_grid - D_comp - h;

[x2_grid,y2_grid,z2_grid] = sph2cart(lon_grid,lat_grid,r2_grid);

lmcosi_gt2_isos=Topo2Grav(flipud(r2_grid'),Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

V1 = Mesh2Volume(x_grid,y_grid,z_grid);
V2 = Mesh2Volume(x2_grid,y2_grid,z2_grid);

M1_isos = rho1_Jh(ind)*V1;
M2_isos = (rho2_Jh(ind)-rho1_Jh(ind))*V2;

w_isos = [M1_isos/(M1_isos+M2_isos) M2_isos/(M1_isos+M2_isos)];
    
lmcosi_gt_isos = WeightSumExpansion(w,{lmcosi_gt1,lmcosi_gt2_isos});

MaxDegreeIsos=min([lmcosi_g(end,1) lmcosi_gt_isos(end,1)]);
lmcosi_isos = lmcosi_g;
lmcosi_isos(:,3:4) = lmcosi_isos(:,3:4) - lmcosi_gt_isos(:,3:4);

% lmcosi_isos(2:6,3:4)=0;

[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_isos,xref,yref,zref);
[gisos_up,gisos_east,gisos_north]=GravityComponents(...
    ax,ay,az,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gisos_up*1e5,'ISOS.dat');

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,gisos_up*1e5); shading interp;
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Isostatic anomaly [mGal]','FontSize',20);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gisos_up*1e5,'ISOS.dat');


[sdl_isos,l_isos] = plm2spec(lmcosi_isos);
[sdl_ba,l_ba] = plm2spec(lmcosi_ba);

figure(h_grav_spec)
h_isos_spec = plot(l_isos,sdl_isos,'om-','LineWidth',3,'MarkerSize',5,'MarkerFaceColor','m');
h_ba_spec = plot(l_ba,sdl_ba,'ob-','LineWidth',3,'MarkerSize',5,'MarkerFaceColor','b');

%% Compute subsurface interface

lmcosi_sub = FindSubRelief(...
    lmcosi_g,lmcosi_t,GM,Rref,rho1_Jh(ind),rho2_Jh(ind),r2_Jh(ind),T);

%% Crustal thickness map

[ri2_sub,lon,lat] = plm2xyz(lmcosi_sub,step);
[ri1,lon,lat]     = plm2xyz(lmcosi_t,step);
[lon,lat] = meshgrid(lon,lat);
ct = (ri1 - ri2_sub)/1000;

image_ct = (ct - min(ct(:)))/...
    (max(ct(:)) - min(ct(:)))*(2^16-1);

imwrite(uint16(image_ct),['Ceres_ct_BW_' ...
    num2str(min(ct(:))) '_' num2str(max(ct(:))) '.png']);

AGUaxes;
pcolorm(lat,lon,ct);
cbar = colorbar('FontSize',20);
ylabel(cbar,'Crustal thickness [km]','FontSize',20);
PrintWhite([fig_folder 'Fig_CT.jpg']);
WriteXYZ(lon,lat,ct,'CT.dat');

%% Anomaly animation

% fig_anim = figure('Position',[1 1 1000 1000]);
% hold on;
% 
% fig_sub_2l  = subplot('Position',[0.5 0.6 0.45 0.35]);
% set(gca, 'FontSize',fntsize);
% hold on;grid on; box on;
% 
% pl_shell     = plot(rho1_Jh,(r1-r2_Jh)/1000,'-k','LineWidth',2);
% pl_shell_pnt = plot(rho1_Jh(1),r1-r2_Jh(1)/1000,'or','MarkerSize',10);
% 
% xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
% ylabel('Shell thickness [km]','FontSize',fntsize);
% 
% fig_sub_map = subplot('Position',[0.05 0.05 0.9 0.45]);
% 
% ax_sub_map=axesm('mollweid','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],...
%     'FontSize',24,'Grid','on','MLabelParallel',...
%     'equator','AngleUnits','degrees','LabelUnits','degrees','ParallelLabel'...
%     , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,...
%     'FontColor',[0 0 0]...
%     ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0],...
%     'FEdgeColor',[0 0 0],'Frame','on');
% 
% axis off;
% 
% fig_sub_2l  = subplot('Position',[0.05 0.6 0.4 0.35]);
% set(gca, 'FontSize',fntsize);
% 
% ax_sub_map2=axesm('mollweid','frame','on','FEdgeColor',[1 1 1],'origin',[0 0],...
%     'FontSize',24,'Grid','on','MLabelParallel',...
%     'equator','AngleUnits','degrees','LabelUnits','degrees','ParallelLabel'...
%     , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',0.5,...
%     'FontColor',[0 0 0]...
%     ,'FontSize',fntsize_sm,'GAltitude',Inf,'Geoid',[1 0],...
%     'FEdgeColor',[0 0 0],'Frame','on');
% 
% writerObj = VideoWriter('BA_Ceres2.avi');
% open(writerObj);
% 
% for i=1:numel(r2_Jh)  
%     
%     [a_Jh(i),~,c_Jh(i)]=fr2abc(r2_Jh(i),fp2_Jh(i),0);
%     lmcosi_gt2=SHRotationalEllipsoid(a_Jh(i),c_Jh(i),MaxDegreeGrav,Rref);
%     w = [M1_Jh(i)/M M2_Jh(i)/M];
%     
%     lmcosi_gt = WeightSumExpansion(w,{lmcosi_gt1,lmcosi_gt2});
%     
%     [ax_gt,ay_gt,az_gt]=GravityAcceleration(...
%         GM,Rref,lmcosi_gt,xref,yref,zref);
%     
%     [g_up_gt,g_east_gt,g_north_gt]=GravityComponents(...
%         ax_gt,ay_gt,az_gt,xref,yref,zref,aref,cref);
%     
% %     WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,(g_up - g_up_gt)*1e5,'BA.dat');
%     
%     axes(ax_sub_map);
%     pcolorm(lat_grid*180/pi,lon_grid*180/pi,(g_up - g_up_gt)*1e5); 
%     caxis([-350 350]);
%     shading interp;
%     cbar = colorbar('FontSize',fntsize);
%     ylabel(cbar,'Bouguer anomaly [mGal]','FontSize',fntsize);
%     drawnow;
%     
%     axes(ax_sub_map2);
% 
%     lmcosi_sub = FindSubRelief(...
%         lmcosi_g,lmcosi_t,GM,Rref,rho1_Jh(i),rho2_Jh(i),r2_Jh(i),T);
%        
%     [ri2_sub,lon,lat] = plm2xyz(lmcosi_sub,step);  
%     [lon,lat] = meshgrid(lon,lat);
%     pcolorm(lat,lon,(ri1 - ri2_sub)/1000);
%     cbar = colorbar('FontSize',fntsize_sm);
%     ylabel(cbar,'Crustal thickness [km]','FontSize',fntsize_sm);
%     
%     title(['$r_{2} = ' num2str(r2_Jh(i)/1000,'%6.2f') ' [km]$; '...
%         '$\rho_{2} = ' num2str(rho2_Jh(i),'%6.2f') ' [kg/m^{3}]$; '...
%         '$\rho_{1} = ' num2str(rho1_Jh(i),'%6.2f') ' [kg/m^{3}]$'],...
%         'FontSize',fntsize,'interpreter','latex');
%     
%     set(pl_shell_pnt,'XData',rho1_Jh(i),'YData',(r1-r2_Jh(i))/1000);
%     
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% 
% close(writerObj);

% AGUaxes;
% pcolorm(lat_grid,lon_grid,g_up_gt); shading interp;
% colorbar('FontSize',fntsize);

%% Admittance

lmcosi_g1_hydro=SHRotationalEllipsoid(aref,cref,MaxDegreeGrav,Rref);   
lmcosi_g_hydro = WeightSumExpansion(w,{lmcosi_g1_hydro,lmcosi_gt2});

r_ell = TriEllRadVec(lat_grid,lon_grid,aref,aref,cref,'rad');
lmcosi_t_hydro = xyz2plm(flipud(r_ell'),MaxDegreeTopo);

lmcosi_t_nonhydro = lmcosi_t;
lmcosi_t_nonhydro(:,3:4) = lmcosi_t_nonhydro(:,3:4) - lmcosi_t_hydro(:,3:4);

lmcosi_g_nonhydro = lmcosi_g;
lmcosi_g_nonhydro(:,3:4) = lmcosi_g_nonhydro(:,3:4) - lmcosi_g_hydro(:,3:4);

lmcosi_gt_nonhydro = lmcosi_gt;
lmcosi_gt_nonhydro(:,3:4) = lmcosi_gt_nonhydro(:,3:4) - lmcosi_g_hydro(:,3:4);

lmcosi_gt_isos_nonhydro = lmcosi_gt_isos;
lmcosi_gt_isos_nonhydro(:,3:4) = lmcosi_gt_isos_nonhydro(:,3:4) - lmcosi_g_hydro(:,3:4);

% [n,Z_gt1] = SphericalHarmonicAdmittance(lmcosi_gt1,lmcosi_t,GM,Rref);

% lmcosi_g_nonhydro(4,3) = 0;
% lmcosi_t_nonhydro(4,3) = 0;
% 
% lmcosi_g(4,3) = 0;
% lmcosi_t(4,3) = 0;

% [n,Z_gw] = SphericalHarmonicAdmittance(lmcosi_g,lmcosi_t,GM,Rref);
[n,Z_gt] = SphericalHarmonicAdmittance(lmcosi_gt_nonhydro,lmcosi_t,GM,Rref);
[n,Z_g] = SphericalHarmonicAdmittance(lmcosi_g_nonhydro,lmcosi_t_nonhydro,GM,Rref);

R0    = lmcosi_t(1,3);
Z_theo  = 3./(2*(n)+1).*((R0/Rref).^n).*(n+1)*GM/(Rref^2)*1e5/R0*1000*...
    rho1_Jh(ind)/rhomean;

[n,Z_isos] = SphericalHarmonicAdmittance(lmcosi_gt_isos_nonhydro,lmcosi_t,GM,Rref);
% [n,Z_isos_lin] = SphericalHarmonicAdmittance(lmcosi_gtisos_lin,lmcosi_t,GM,Rref);

fig_Z=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XTick',1:100);

plot(n,Z_g,'-b','LineWidth',3,'MarkerSize',5);
plot(n,Z_gt,'-k','LineWidth',3,'MarkerSize',5);
plot(n,Z_isos,'-g','LineWidth',3,'MarkerSize',5);

save('Z_g_obs.mat','n','Z_g');

% plot(n,Z_theo,'-r','LineWidth',3,'MarkerSize',5);
% plot(n,Z_isos,'-g','LineWidth',3,'MarkerSize',5);
% plot(n,Z_isos_lin,'-m','LineWidth',3,'MarkerSize',5);

legend({'Observed','2-layer','Isostatic'},'FontSize',fntsize_sm);
xlabel('Degree','FontSize',fntsize);
ylabel('Admittance [mGal/km]','FontSize',fntsize);

PrintWhite(fig_Z,[fig_folder 'Fig_Z.jpg']);

%% Correlation
cor_isos = SphericalHarmonicCorrelation(lmcosi_gt_isos_nonhydro,lmcosi_gt_nonhydro);
cor_gt = SphericalHarmonicCorrelation(lmcosi_gt_nonhydro,lmcosi_gt_nonhydro);
cor_g = SphericalHarmonicCorrelation(lmcosi_g_nonhydro,lmcosi_gt_nonhydro);

fig_R=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XTick',1:100);

xlim([2 MaxDegreeGrav]);
% ylim([-1 1]);

n = 0:MaxDegreeGrav;
h_r_g = plot(n,cor_g,'-b','LineWidth',3,'MarkerSize',5);
h_r_gt = plot(n,cor_gt,'-r','LineWidth',3,'MarkerSize',5);
h_r_gt_isos = plot(n,cor_isos,'-g','LineWidth',3,'MarkerSize',5);

legend([h_r_g h_r_gt h_r_gt_isos], {'Observed','2-layer','Isostatic'},'FontSize',fntsize_sm);
xlabel('Degree','FontSize',fntsize);
ylabel('Correlation','FontSize',fntsize);

PrintWhite(fig_R,[fig_folder 'Fig_R.jpg']);

%% Effective density

[sdl_g,l_g] = plm2spec(lmcosi_g);
[sdl_gt,l_gt] = plm2spec(lmcosi_gt);

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

plot(l_g(l_g>1),sdl_g(l_g>1)./sdl_gt(l_g>1)*rhomean,'-k','LineWidth',3,'MarkerSize',5);

xlabel('Degree','FontSize',fntsize);
ylabel('Effective density [kg/m^3]','FontSize',fntsize);

