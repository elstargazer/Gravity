ccc

%
G = 6.67408e-11;
GM = 62.6315e9;
%
fntsize = 12;
fntsize_sm = 10;
im_size = [0 0 13 9];

%%
step = 0.1;
Npts = 100000;

%% Load SPC shape

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150828_512.bds';
fig_folder='~/Dawn/Papers/CeresPaper1/';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

[x_grid_spc,y_grid_spc,z_grid_spc]=ReadSPC(full_filename,step,'grid');
[x_rand_spc,y_rand_spc,z_rand_spc]=ReadSPC(full_filename,Npts,'rand');

r_grid_spc=sqrt(x_grid_spc.^2+y_grid_spc.^2+z_grid_spc.^2);

r_rand_vec_spc = [x_rand_spc' y_rand_spc' z_rand_spc'];
r_grid_vec_spc = [x_grid_spc(:) y_grid_spc(:) z_grid_spc(:)];

%% Load SPG model

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPG/Survey/';
shape_filename='global.bds';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

[x_grid_spg,y_grid_spg,z_grid_spg]=ReadSPC(full_filename,step,'grid');
[x_rand_spg,y_rand_spg,z_rand_spg]=ReadSPC(full_filename,Npts,'rand');

x_grid_spg=x_grid_spg/1000;
y_grid_spg=y_grid_spg/1000;
z_grid_spg=z_grid_spg/1000;

x_rand_spg=x_rand_spg/1000;
y_rand_spg=y_rand_spg/1000;
z_rand_spg=z_rand_spg/1000;

r_grid_spg=sqrt(x_grid_spg.^2+y_grid_spg.^2+z_grid_spg.^2);

r_rand_vec_spg = [x_rand_spg' y_rand_spg' z_rand_spg'];
r_grid_vec_spg = [x_grid_spg(:) y_grid_spg(:) z_grid_spg(:)];

[B,L,H] = XYZ2BLH(x_grid_spc,y_grid_spc,z_grid_spc,482.000,Eccentricity(482.000,446.000));

%%

[lon_grid,lat_grid,r_grid]=cart2sph(x_grid_spc,y_grid_spc,z_grid_spc);
[lon_rand,lat_rand,r_rand]=cart2sph(x_rand_spc,y_rand_spc,z_rand_spc);


%% Compute the best-fit ellipsoid

[ center_spc, ell_spc, evecs_spc, v_spc ] = ellipsoid_fit( r_rand_vec_spc, 0 );
[ center_spg, ell_spg, evecs_spg, v_spg ] = ellipsoid_fit( r_rand_vec_spg, 0 );

[lambda_ax_spc,fi_ax_spc]=AxesOrientation(evecs_spc);
[lambda_ax_spg,fi_ax_spg]=AxesOrientation(evecs_spg);

H_grid_spc = EllipsoidHeight([x_grid_spc(:) y_grid_spc(:) z_grid_spc(:)],...
    center_spc,ell_spc,evecs_spc);
H_rand_spc = EllipsoidHeight([x_rand_spc(:) y_rand_spc(:) z_rand_spc(:)],...
    center_spc,ell_spc,evecs_spc);
H_rand_rms_spc = sqrt(sum(H_rand_spc.*H_rand_spc)/numel(H_rand_spc));

H_grid_spc=reshape(H_grid_spc,size(x_grid_spc));

H_grid_spg = EllipsoidHeight([x_grid_spg(:) y_grid_spg(:) z_grid_spg(:)],...
    center_spg,ell_spg,evecs_spg);
H_rand_spg = EllipsoidHeight([x_rand_spg(:) y_rand_spg(:) z_rand_spg(:)],...
    center_spg,ell_spg,evecs_spg);
H_rand_rms_spg = sqrt(sum(H_rand_spg.*H_rand_spg)/numel(H_rand_spg));

H_grid_spg=reshape(H_grid_spg,size(x_grid_spg));

%% 

AGUaxes;
title('SPC');
pcolorm(lat_grid*180/pi,lon_grid*180/pi,H_grid_spc); 
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Elevation [km]','FontSize',fntsize);

AGUaxes;
title('SPG');
pcolorm(lat_grid*180/pi,lon_grid*180/pi,H_grid_spg); 
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Elevation [km]','FontSize',fntsize);

%% 

AGUaxes;
title('SPC - SPG');
pcolorm(lat_grid*180/pi,lon_grid*180/pi,H_grid_spc-H_grid_spg); 
cbar = colorbar('FontSize',fntsize);
ylabel(cbar,'Elevation difference [km]','FontSize',fntsize);

%% Elevation Histogram

figure;
hold on;box on;grid on;
set(gcf, 'Units','centimeters', 'Position',im_size);
set(gcf, 'PaperPositionMode','auto');
set(gca, 'FontSize',fntsize);

hist(H_rand_spg,20,'b');
hist(H_rand_spc,20,'r');

alpha(0.5);

xlabel('Elevation [km]','FontSize',fntsize);
ylabel('Number of points []','FontSize',fntsize);

% 
figure;
hold on;box on;grid on;
set(gcf, 'Units','centimeters', 'Position',im_size);
set(gcf, 'PaperPositionMode','auto');
set(gca, 'FontSize',fntsize);

hist(H_grid_spc(:)-H_grid_spg(:),20,'k');
xlabel('Elevation [km]','FontSize',fntsize);
ylabel('Number of points []','FontSize',fntsize);
PrintWhite([fig_folder 'Fig_hist.eps']);

%% Spherical harmonics

L=540;

lmcosi_ceres_spc = xyz2plm(r_grid_spc',L);
lmcosi_ceres_spg = xyz2plm(r_grid_spg',L);

%% Power spectra

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;
set(gca,'XScale','log');
set(gca,'YScale','log');

lmcosi_ceres_diff = lmcosi_ceres_spg;
lmcosi_ceres_diff(:,3:4) = lmcosi_ceres_spc(:,3:4) - lmcosi_ceres_spg(:,3:4);
[sdl_spc,l_spc] = plm2spec(lmcosi_ceres_spc);
[sdl_spg,l_spg] = plm2spec(lmcosi_ceres_spg);
[sdl_diff,l_diff] = plm2spec(lmcosi_ceres_diff);

plot(l_spc,sdl_spc,'r-');
plot(l_spg,sdl_spg,'b-');
plot(l_diff,sdl_diff,'k-');

legend({'SPC','SPG','SPC - SPG'},'FontSize',fntsize_sm);

xlabel('Degree','FontSize',fntsize);
ylabel('Topography power','FontSize',fntsize);

%% Correlation

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;
corr=SphericalHarmonicCorrelation(lmcosi_ceres_spc,lmcosi_ceres_spg);

plot(corr,'-k')

xlabel('Spherical harmonics degree','FontSize',fntsize);
ylabel('Correlation','FontSize',fntsize);

PrintWhite([fig_folder 'Fig_Corr.jpg']);


%% Isotropy


figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;
q_spc=IsotropicRatio(lmcosi_ceres_spc,lmcosi_ceres_spc);
q_spg=IsotropicRatio(lmcosi_ceres_spg,lmcosi_ceres_spg);

plot(q_spc,'-r')
plot(q_spg,'-b')

ylim([0 3]);

legend({'SPC','SPG'},'FontSize',fntsize_sm);


xlabel('Spherical harmonics degree','FontSize',fntsize);
ylabel('Isotropic ratio','FontSize',fntsize);

PrintWhite([fig_folder 'Fig_IsotropicRatio.jpg']);


%% Compare gravity 





