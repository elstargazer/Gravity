ccc

fntsize = 12;
fntsize_sm = 10;
im_size = [0 0 13 9];

%% Load OpNav shape

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150828_512.bds';
fig_folder='~/Dawn/Papers/CeresPaper1/';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

step = 0.1;
Npts = 100000;
[x_grid,y_grid,z_grid]=ReadSPC(full_filename,step,'grid');
[x_rand,y_rand,z_rand]=ReadSPC(full_filename,Npts,'rand');

r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);

V = Mesh2Volume(x_grid,y_grid,z_grid);

[lon_grid,lat_grid,r_grid]=cart2sph(x_grid,y_grid,z_grid);
[lon_rand,lat_rand,r_rand]=cart2sph(x_rand,y_rand,z_rand);

cond_lat = (lat_rand > - 90/180*pi);
r_rand_vec = [x_rand(cond_lat)' y_rand(cond_lat)' z_rand(cond_lat)'];
r_grid_vec = [x_grid(:) y_grid(:) z_grid(:)];
%% Compute the best-fit ellipsoid

[ center, ell, evecs, v ] = ellipsoid_fit( r_rand_vec, 0 );
[ center_xy, ell_xy, evecs_xy, v_xy ] = ellipsoid_fit( r_rand_vec, 2 , 'xy' );
[ center_bi, ell_bi, evecs_bi, v_bi ] = ellipsoid_fit( r_rand_vec, 4);

% [ center, radii, evecs, v ] = ellipsoid_fit_err( r_rand, 4 );
% radii_std=std(radii,0,2);
% center_std=std(center,0,2);

[v2,e,d] = ellipsoidfit(x_rand,y_rand,z_rand);
v3 = ellipsoidfit_direct(x_rand,y_rand,z_rand);
v4 = ellipsoidfit_leastsquares(x_rand,y_rand,z_rand);

[center2,ell2,quat2,R2] = ellipsoid_im2ex(v2);
[center3,ell3,quat3,R3] = ellipsoid_im2ex(v3);
[center4,ell4,quat4,R4] = ellipsoid_im2ex(v4);

offset=norm(center);
[lambda_ax,fi_ax]=AxesOrientation(evecs);

[fp,fq,fabp] = ell2f(ell);

Rvol=prod(ell).^(1/3);

H_grid = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center,ell,evecs);

H_rand = EllipsoidHeight([x_rand(:) y_rand(:) z_rand(:)],...
    center,ell,evecs);

H_rand_rms = sqrt(sum(H_rand.*H_rand)/numel(H_rand));

H_grid_xy = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center_xy,ell_xy,evecs_xy);

H_grid_bi = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center_bi,ell_bi,evecs_bi);

H_grid=reshape(H_grid,size(x_grid));
H_grid_xy=reshape(H_grid_xy,size(x_grid));
H_grid_bi=reshape(H_grid_bi,size(x_grid));

[lambda,fi,rad]=cart2sph(x_grid,y_grid,z_grid);

%% Write Ellipsoidal elevation grid

WriteXYZ(lambda*180/pi,fi*180/pi,H_grid,[shape_folder shapename '_gen.xyz']);
WriteXYZ(lambda*180/pi,fi*180/pi,H_grid_xy,[shape_folder shapename '_xy.xyz']);
WriteXYZ(lambda*180/pi,fi*180/pi,r_grid,[shape_folder shapename '_rad.xyz']);
WriteXYZ(lambda*180/pi,fi*180/pi,H_grid_bi,[shape_folder shapename '_bi.xyz']);

%% Elevation Histogram

figure;
hold on;box on;grid on;
set(gcf, 'Units','centimeters', 'Position',im_size);
set(gcf, 'PaperPositionMode','auto');
set(gca, 'FontSize',fntsize);

hist(H_rand,20,'b');

xlabel('Elevation [km]','FontSize',fntsize);
ylabel('Number of points []','FontSize',fntsize);

PrintWhite([fig_folder 'Fig_hist.eps']);

%% Show Ceres

% im = imread('~/Dawn/Reflectance/Ceres/global.png','png');
% ShowCeres

%% Spherical harmonics
% 
% VestaShapeModel='~/Dawn/MATLAB/Vesta20140513shape_geoc_elev_3m.grd';
% [~,~,r_vesta]=ReadGRD(VestaShapeModel);

L=550;

aref=482.788081;
bref=aref;
cref=447.8000;

aref_vesta=281;
bref_vesta=aref_vesta;
cref_vesta=226;

r_ell=TriEllRadVec(fi,lambda,aref,bref,cref,'rad');
% r_ell_vesta=TriEllRadVec(fi,lambda,aref_vesta,bref_vesta,cref_vesta,'rad');

lmcosi_ceres = xyz2plm((r_grid'),L,'im');

%% Produce filtered topography
% lmcosi_ceres_filtered = ApplyCosFilter(lmcosi_ceres,15,30);
% [r_sh,lon_sh,lat_sh] = plm2xyz(lmcosi_ceres_filtered,step);
% [lon_sh,lat_sh] = meshgrid(lon_sh,lat_sh);
% [x_sh,y_sh,z_sh] = sph2cart(lon_sh/180*pi,lat_sh/180*pi,flipud(r_sh));
% H_grid_sh = EllipsoidHeight([x_sh(:) y_sh(:) z_sh(:)],...
%     center_xy,ell_xy,evecs_xy);
% 
% WriteXYZ(lon_sh,lat_sh,H_grid_sh,[shape_folder shapename '_gen_sh.xyz']);
% 
% H_grid_sh_r = reshape(H_grid_sh,size(x_sh));
% H_grid_r = reshape(H_grid_xy,size(x_grid));
% 
% image_sh =  (H_grid_sh_r - min(H_grid_sh_r(:)))/...
%     (max(H_grid_sh_r(:)) - min(H_grid_sh_r(:)))*(2^8-1);
% 
% image_grid = (H_grid_r - min(H_grid_r(:)))/...
%     (max(H_grid_r(:)) - min(H_grid_r(:)))*(2^8-1);
% map = jet(256);
% 
% imwrite(uint8(image_sh),map,'CeresFilteredColor_xy2.png');
% imwrite(uint8(flipud(image_grid')),map,'Ceres_Color_xy.png');
% 
% 
% figure; hold on; colormap jet
% pcolor(flipud(H_grid_sh_r)); shading interp;
% cbar = colorbar('FontSize',20);
% ylabel(cbar,'Elevation [km]','FontSize',20);
% 
% 
% figure; hold on; colormap jet
% pcolor(flipud(H_grid_r')); shading interp;
% cbar = colorbar('FontSize',20);
% ylabel(cbar,'Elevation [km]','FontSize',20);
% 
% 
% image_sh = (H_grid_sh_r - min(H_grid_r(:)))/...
%     (max(H_grid_r(:)) - min(H_grid_r(:)))*(2^16-1);
% 
% image_grid = (H_grid_r - min(H_grid_r(:)))/...
%     (max(H_grid_r(:)) - min(H_grid_r(:)))*(2^16-1);
% 
% imwrite(uint16(image_sh),'CeresFilteredBW_xy.png');
% imwrite(uint16(flipud(image_grid')),'CeresBW_xy.png');
% 
% [min(H_grid_sh_r(:)) max(H_grid_sh_r(:))]
% [min(H_grid_r(:)) max(H_grid_r(:))]
% 
% % lmcosi_ell = xyz2plm(r_ell',L);
% % lmcosi_ell_vesta = xyz2plm(r_ell_vesta',L);
% % lmcosi_vesta = xyz2plm(r_vesta,L);
% 
% % lmcosi(:,3:4)=lmcosi(:,3:4)/lmcosi(1,3);
% % lmcosi_ell(:,3:4)=lmcosi_ell(:,3:4)/lmcosi_ell(1,3);
% % lmcosi_vesta(:,3:4)=lmcosi_vesta(:,3:4)/lmcosi_vesta(1,3);
% % lmcosi_ell_vesta(:,3:4)=lmcosi_ell_vesta(:,3:4)/lmcosi_ell_vesta(1,3);

%% Spectra

[sdl,l]=plm2spec(lmcosi_ceres);
[r_grid_sh,lat_sh,lon_sh] = plm2xyz(lmcosi_ceres,step);

figure; hold on;
set(gca,'FontSize',fntsize);
set(gca,'YScale','log');

plot(l,sdl,'ro-','LineWidth',2,'MarkerSize',5,...
    'MarkerFaceColor','r');

% plot(l_ell(1:2:end),sdl_ell(1:2:end),...
%     'r*--','LineWidth',2,'MarkerSize',5);

% plot(l_vesta,sdl_vesta,'bo-','LineWidth',3,'MarkerSize',20,'MarkerFaceColor','b');
% plot(l_ell_vesta(1:2:end),sdl_ell_vesta(1:2:end),'b*--','LineWidth',3,'MarkerSize',20);

% legend({'Ceres shape','Ceres hydrostatic shape','Vesta','Vesta hydrostatic shape'},'FontSize',20);

xlabel('Degree []','FontSize',fntsize);
ylabel('Spectral density  []','FontSize',fntsize);
set(gcf, 'Units','centimeters', 'Position',im_size)
box on;

% xlim([0 10]);

% q=IsotropicRatio(lmcosi(2:end,:),lmcosi(2:end,:));
% 
% figure; hold on;
% plot(1:numel(q),q,'-bo')
% set(gca,'Yscale','log');


%%
lmcosi_ceres(:,3:4)=lmcosi_ceres(:,3:4)*1000;
lmcosi_ceres(4,3)=0;
[k_Ceres,sdl_Ceres]=PowerSpectrum(lmcosi_ceres);

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

set(gca,'FontSize',fntsize);
set(gca,'XScale','log');
set(gca,'YScale','log');

h = plot(k_Ceres,sdl_Ceres,'-k');

[pfit,S] = polyfit(log10(k_Ceres),...
    log10(sdl_Ceres),1);
[log_sdl_Ceres_err,log_delta_sdl] = polyconf(pfit,log10(k_Ceres),S);

h_fit1 = plot(k_Ceres,10.^log_sdl_Ceres_err,'-r','LineWidth',1);

plot(k_Ceres,10.^(log_sdl_Ceres_err+log_delta_sdl),...
    '--r','LineWidth',1);
plot(k_Ceres,10.^(log_sdl_Ceres_err-log_delta_sdl),...
    '--r','LineWidth',1);

ind_fit=(k_Ceres>4e-3);

[pfit,S] = polyfit(log10(k_Ceres(ind_fit)),...
    log10(sdl_Ceres(ind_fit)),1);
[log_sdl_Ceres_err,log_delta_sdl] = polyconf(pfit,log10(k_Ceres),S);

% h_fit2 = plot(k_Ceres,10.^log_sdl_Ceres_err,'-b','LineWidth',1);
% 
% plot(k_Ceres,10.^(log_sdl_Ceres_err+log_delta_sdl),...
%     '--b','LineWidth',1);
% plot(k_Ceres,10.^(log_sdl_Ceres_err-log_delta_sdl),...
%     '--b','LineWidth',1);

xlabel('Frequency [cycles/km]','FontSize',fntsize);
ylabel('Topography power [km^3]','FontSize',fntsize);

% xlim([5*1e-1 2e2]);
set(gca,'XTick',10.^(-5:1:0));
set(gca,'YTick',10.^(0:2:16));

xlim([1e-4 1e0]);
ylim([1e0 1e15]);
% legend([h h2],{'Ceres JPL','Ceres DLR'},'FontSize',fntsize_sm);

legend({'Ceres spectrum','Power law'},'FontSize',fntsize_sm);
PrintWhite([fig_folder 'Fig_spectrum_fits.jpg']);






