ccc

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

%% Load OpNav shape

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150702_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150702_256.bds';
fig_folder='~/Dawn/Papers/CeresPaper1/';

[~,shapename,~] = fileparts(shape_filename) ;

full_filename = [shape_folder shape_filename];

step = 0.1;
Npts = 100000;
[x_grid,y_grid,z_grid]=ReadSPC(full_filename,step,'grid');
[x_rand,y_rand,z_rand]=ReadSPC(full_filename,Npts,'rand');

r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);
r_rand_vec = [x_rand' y_rand' z_rand'];
V = Mesh2Volume(x_grid,y_grid,z_grid);

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

L=360;

aref=482.788081;
bref=aref;
cref=447.8000;

aref_vesta=281;
bref_vesta=aref_vesta;
cref_vesta=226;

r_ell=TriEllRadVec(fi,lambda,aref,bref,cref,'rad');
% r_ell_vesta=TriEllRadVec(fi,lambda,aref_vesta,bref_vesta,cref_vesta,'rad');

lmcosi_ceres = xyz2plm(r_grid',L,'im');
lmcosi_ell = xyz2plm(r_ell',L);

% lmcosi_ell_vesta = xyz2plm(r_ell_vesta',L);
% lmcosi_vesta = xyz2plm(r_vesta,L);


% lmcosi(:,3:4)=lmcosi(:,3:4)/lmcosi(1,3);
% lmcosi_ell(:,3:4)=lmcosi_ell(:,3:4)/lmcosi_ell(1,3);
% lmcosi_vesta(:,3:4)=lmcosi_vesta(:,3:4)/lmcosi_vesta(1,3);
% lmcosi_ell_vesta(:,3:4)=lmcosi_ell_vesta(:,3:4)/lmcosi_ell_vesta(1,3);

[sdl,l]=plm2spec(lmcosi_ceres);
[sdl_ell,l_ell]=plm2spec(lmcosi_ell);
% [sdl_vesta,l_vesta]=plm2spec(lmcosi_vesta);
% [sdl_ell_vesta,l_ell_vesta]=plm2spec(lmcosi_ell_vesta);

figure; hold on;
set(gca,'FontSize',fntsize);
set(gca,'YScale','log');

plot(l,sdl,'ro-','LineWidth',2,'MarkerSize',5,...
    'MarkerFaceColor','r');

plot(l_ell(1:2:end),sdl_ell(1:2:end),...
    'r*--','LineWidth',2,'MarkerSize',5);
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

% [k_Ceres,sdl_Ceres]=PowerSpectrum(lmcosi_ceres);
lmcosi_ceres(:,3:4)=lmcosi_ceres(:,3:4)*1000;
[k_Ceres,sdl_Ceres]=PowerSpectrum(lmcosi_ceres);

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

set(gca,'FontSize',fntsize);
set(gca,'XScale','log');
set(gca,'YScale','log');

plot(k_Ceres,sdl_Ceres,'-k');

xlabel('Frequency [cycles/km]','FontSize',fntsize);
ylabel('Topography power [km^3]','FontSize',fntsize);

% xlim([5*1e-1 2e2]);
set(gca,'XTick',10.^(-5:1:0));
set(gca,'YTick',10.^(0:2:16));

xlim([1e-5 1e0]);
ylim([1e0 1e15]);

legend({'Ceres'},'FontSize',fntsize_sm);

PrintWhite([fig_folder 'Fig_spectrum.jpg']);












