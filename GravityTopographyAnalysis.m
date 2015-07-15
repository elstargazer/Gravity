ccc

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

%% input parameters
shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150702_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150702_256.bds';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

GM = 62.6253e9;
G = 6.67384e-11;
Rref=476000;
aref=481000;
cref=446000;
step = 0.5;
r1 = 470000;
T = 9.073859324514187; % DLR
Npts = 50;

rhomean=2150;
M=GM/G;

MaxDegreeTopo = 100;
MaxDegreeGrav = 2;
MaxTopoPower  = 4;

%% get shape model in SH

[x_grid,y_grid,z_grid]=ReadSPC(full_filename,step,'grid');

x_grid=x_grid*1000;
y_grid=y_grid*1000;
z_grid=z_grid*1000;

[lon_grid,lat_grid,r_grid]=cart2sph(x_grid,y_grid,z_grid);
% reference ellipsoid surface
[xref,yref,zref]=TriEllRadVec(lat_grid,lon_grid,aref,aref,cref,'xyz');

eccref=Eccentricity(aref,cref);
[B,L,H]=XYZ2BLH(x_grid,y_grid,z_grid,aref,eccref);

% lmcosi_topo = xyz2plm(r_grid,L);

%% get gravity model in SH

lmcosi_g = [0 0 1 0;
    1 0 0 0;
    1 1 0 0;
    2 0 -1.14e-2 0;
    2 1 -8.65e-7 -1.43e-5;
    2 2 2.34e-4 -2.71e-4];

   lmcosi_gt_ryan = [0 0 1 0;
    1 0 0 0;
    1 1 0 0;
    2 0 -1.14e-29 0;
    2 1 -8.65e-79 -1.43e-59;
    2 2 -0.291084874218D-04 0.716547885022D-03];


     

J2obs = -lmcosi_g(4,3);

%% plot topography

AGUaxes;
pcolorm(lat_grid,lon_grid,H);

%% hydrostatic gravity

% Grid of core radii and densities
% rcore=2500:500:470000;
% rhocoreg=2000:100:4000;

r2=linspace(10000,470000,Npts);
rho2=linspace(rhomean,5000,Npts);

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

r2_Jh   = CJhyd(1,2:end);
rho2_Jh = CJhyd(2,2:end);

rho1_Jh = griddata(r2i,rho2i,rho1i,r2_Jh,rho2_Jh,'linear');
M2_Jh   = griddata(r2i,rho2i,M2,r2_Jh,rho2_Jh,'linear');
f2_Jh   = griddata(r2i,rho2i,f2i,r2_Jh,rho2_Jh,'linear');

M1_Jh = M - M2_Jh;

fig1=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on;

plot(rho1_Jh,(r1-r2_Jh)/1000)

% we have r1, r2_Jh, rho1_Jh, rho2_Jh = family of solutions for J2
%% plot gravity

% compute just gravity 
[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_g,xref,yref,zref);
[g_up,g_east,g_north]=GravityComponents(ax,ay,az,xref,yref,zref,aref,cref);

% computing free-air anomaly
lmcosi_fa = lmcosi_g;
lmcosi_fa(4,3) = 0;
lmcosi_fa(1,3) = 0;

[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_fa,xref,yref,zref);
[g_up_fa,g_east_fa,g_north_fa]=GravityComponents(ax,ay,az,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,g_up_fa*1e5,'FA.dat');

AGUaxes;
pcolorm(lat_grid,lon_grid,g_up_fa*1e5); shading interp;
colorbar('FontSize',fntsize);

%% compute gravity from shape

% lmcosi_gt1=TopoSH2GravitySH(flipud(r_grid'),GM,rhomean,Rref,...
%     MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_gt1=Topo2Grav(flipud(r_grid'),Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

%% model gravity
% compute gravity from shape
[ax_gt1,ay_gt1,az_gt1]=GravityAcceleration(...
    GM,Rref,lmcosi_gt1,xref,yref,zref);
[g_up_gt1,g_east_gt1,g_north_gt1]=GravityComponents(...
    ax_gt1,ay_gt1,az_gt1,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,g_up_gt1*1e5,'GT.dat');

% compute gravity from core

a_Jh = zeros(size(r2_Jh));
c_Jh = a_Jh;

AGUaxes;

for i=1:numel(r2_Jh)   
    [a_Jh(i),~,c_Jh(i)]=fr2abc(r2_Jh(i),f2_Jh(i),0);
    lmcosi_gt2=SHRotationalEllipsoid(a_Jh(i),c_Jh(i),MaxDegreeGrav,Rref); 
    w = [M1_Jh(i)/M M2_Jh(i)/M];
    
    lmcosi_gt = WeightSumExpansion(w,{lmcosi_gt1,lmcosi_gt2});  
    
    [ax_gt,ay_gt,az_gt]=GravityAcceleration(...
        GM,Rref,lmcosi_gt,xref,yref,zref);
    
    [g_up_gt,g_east_gt,g_north_gt]=GravityComponents(...
        ax_gt,ay_gt,az_gt,xref,yref,zref,aref,cref);
    
    WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,(g_up - g_up_gt)*1e5,'BA.dat');
       
    pcolorm(lat_grid*180/pi,lon_grid*180/pi,(g_up - g_up_gt)*1e5); shading interp;
    cbar = colorbar('FontSize',fntsize);
    ylabel(cbar,'Bouguer anomaly [mGal]','FontSize',fntsize);
    drawnow;
    pause(0.1); 
    
end

% AGUaxes;
% pcolorm(lat_grid,lon_grid,g_up_gt); shading interp;
% colorbar('FontSize',fntsize);

%% compute Bouguer anomaly

MaxDegreeBouguer=min([lmcosi_g(end,1) lmcosi_gt(end,1)]);
lmcosi_ba = lmcosi_g;
lmcosi_ba(:,3:4) = lmcosi_ba(:,3:4) - lmcosi_gt(:,3:4);
lmcosi_ba(4,3)=0;

[ax,ay,az]=GravityAcceleration(GM,Rref,lmcosi_ba,xref,yref,zref);
[gba_up,gba_east,gba_north]=GravityComponents(...
    ax,ay,az,xref,yref,zref,aref,cref);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gba_up*1e5,'BA.dat');

AGUaxes;
pcolorm(lat_grid,lon_grid,gba_up*1e5); shading interp;
colorbar('FontSize',fntsize);


lmcosi_g_noJ2 = lmcosi_g;
lmcosi_gt_noJ2 = lmcosi_gt;
lmcosi_g_noJ2(4,3)=0;
lmcosi_gt_noJ2(4,3)=0;

cor = SphericalHarmonicCorrelation(lmcosi_g,lmcosi_gt);
cor_noJ2 = SphericalHarmonicCorrelation(lmcosi_g_noJ2,lmcosi_gt_noJ2);
cor_noJ2_2 = SphericalHarmonicCorrelation(lmcosi_g_noJ2,lmcosi_gt_ryan);


% Z_noJ2 = SphericalHarmonicAdmittance(lmcosi_g_noJ2,lmcosi_gt_noJ2);


%% isostatic anomaly


%% Compute sub-surface relief
