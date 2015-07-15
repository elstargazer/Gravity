ccc

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

%% Input parameters

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150702_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150702_256.bds';
step = 0.5;
Rref = 476000;
D = 30000;
rho_crust = 1000;
rho_mantle = 2300;
nmaxt = 100;
nmaxgt = 20;
hmaxt = 1;
GM = 62.6253e9;
G = 6.67384e-11;

aref=481000;
cref=446000;

%% Read shape data

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

[x_grid,y_grid,z_grid]=ReadSPC(full_filename,step,'grid');
[lon_grid,lat_grid,r_grid] = cart2sph(x_grid,y_grid,z_grid);
r_grid=flipud(r_grid')*1000;

[xref,yref,zref]=TriEllRadVec(lat_grid,lon_grid,aref,aref,cref,'xyz');

V = Mesh2Volume(x_grid,y_grid,z_grid);
rho_mean = GM/G/V/1e9;

%% Compute gravity

lmcosi_shape = xyz2plm(r_grid,nmaxt,'im');
R0 = lmcosi_shape(1,3);

lmcosi_gravtopo2 = TopoSH2GravitySH(...
    r_grid,GM,rho_mean,Rref,nmaxt,nmaxgt,hmaxt);

lmcosi_gravtopo = Topo2Grav(...
    r_grid,Rref,nmaxt,nmaxgt,hmaxt);

lmcosi_isosgrav = Topo2IsosGrav(...
    r_grid,Rref,D,rho_crust,rho_mantle,rho_mean,nmaxt,nmaxgt,hmaxt);


[ax_gt,ay_gt,az_gt]=GravityAcceleration(...
    GM,Rref,lmcosi_gravtopo,xref,yref,zref);
[g_up_gt,g_east_gt,g_north_gt]=GravityComponents(...
    ax_gt,ay_gt,az_gt,xref,yref,zref,aref,cref);

AGUaxes;
pcolorm(lat_grid*180/pi,lon_grid*180/pi,g_up_gt);

%% Compute admittance

[n_homo,Z_homo]   = SphericalHarmonicAdmittance(lmcosi_gravtopo,lmcosi_shape,GM,Rref);
[n_homo2,Z_homo2] = SphericalHarmonicAdmittance(lmcosi_gravtopo2,lmcosi_shape,GM,Rref);
[n_isos,Z_isos]   = SphericalHarmonicAdmittance(lmcosi_isosgrav,lmcosi_shape,GM,Rref);

[R_homo]   = SphericalHarmonicCorrelation(lmcosi_gravtopo2,lmcosi_shape);

Z_theo  = 3./(2*(n_homo)+1).*((R0/Rref).^n_homo).*(n_homo+1)*GM/(Rref^2)*1e5/1000;

%% Plot admittances 

fig_iso = figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

h_homo2 = plot(n_homo2,Z_homo2,'-g');
h_homo = plot(n_homo,Z_homo,'--r');
h_isos = plot(n_isos,Z_isos,'-b');
h_theo = plot(n_homo,Z_theo,'-k');

legend([h_homo2 h_homo h_isos h_theo],{'Homo','Homo2','Isostatic','Theor'},...
    'FontSize',fntsize_sm);

xlabel('Spherical harmonic degree','FontSize',fntsize);
ylabel('Admittance []','FontSize',fntsize);