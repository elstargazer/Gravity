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
filename_grav = '/Users/antonermakov/Dawn/CeresGravityModel/CERES08A01/JGC08A01.sha';

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
step  = 0.5;
r1    = 470000;
T     = 9.073859324514187; % DLR
Npts  = 50;

MaxDegreeTopo = 70;
MaxDegreeGrav = 5;
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

lmcosi_t = xyz2plm(flipud(r_grid'),MaxDegreeTopo);

%% get gravity model in SH

[lmcosi_g,Rref,GM,GM_std]=ReadGRAILGravityModel(filename_grav);
lmcosi_g = [0 0 1 0 0 0; lmcosi_g];

lmcosi_g = TruncateGravityModel(lmcosi_g,MaxDegreeGrav,1);
% GM = GM*1e9;
M=GM/G;
rhomean = M/V;

J2obs = -lmcosi_g(4,3);

%% Plot topography with respect to equipotential surface

H_eq = Height2Equipotential(full_filename,lat_grid,lon_grid,GM,Rref,lmcosi_g,T);

%% hydrostatic gravity

r2=linspace(10000,470000,Npts);
r2 = r2(1:end-1);
r2_add = linspace(r2(end),r1,5);
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

rho1_lin = 800:50:2000;
st_lin = interp1(rho1_Jh(~isnan(rho1_Jh)),(r1-r2_Jh(~isnan(rho1_Jh)))/1000,rho1_lin,'cubic');
rho2_lin = interp1(rho1_Jh(~isnan(rho1_Jh)),rho2_Jh(~isnan(rho1_Jh)),rho1_lin,'cubic');
r2_lin = interp1(rho1_Jh(~isnan(rho1_Jh)),r2_Jh(~isnan(rho1_Jh)),rho1_lin,'cubic');

xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);

gi = ginput(1);
ind = find(abs(rho1_Jh - gi(1)) == min(abs(rho1_Jh - gi(1))));
plot(rho1_Jh(ind),(r1-r2_Jh(ind))/1000,'or','MarkerSize',10);

% we have r1, r2_Jh, rho1_Jh, rho2_Jh = family of solutions for J2

%% compute gravity from shape

lmcosi_gt1=Topo2Grav(flipud(r_grid'),Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

%% model gravity
% compute gravity from shape
[ax_gt1,ay_gt1,az_gt1]=GravityAcceleration(...
    GM,Rref,lmcosi_gt1,xref,yref,zref);
[g_up_gt1,g_east_gt1,g_north_gt1]=GravityComponents(...
    ax_gt1,ay_gt1,az_gt1,xref,yref,zref,aref,cref);

a_Jh = zeros(size(r2_Jh));
c_Jh = a_Jh;

%% Compute isostatic anomaly

D_comp = r1 - r2_Jh(ind);
 
t = FindCrustalRoot(r1,D_comp,H_eq,rho1_Jh(ind),rho2_Jh(ind)-rho1_Jh(ind));

r2_grid = TriEllRadVec(lat_grid,lon_grid,a_Jh,a_Jh,c_Jh,'rad');
r2_grid = r2_grid - t;

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

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,gisos_up*1e5); shading interp;
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Isostatic anomaly [mGal]','FontSize',20);

WriteXYZ(lon_grid*180/pi,lat_grid*180/pi,gisos_up*1e5,'ISOS.dat');
