ccc

fntsize=35;
fntsize_sm=23;
imsize=[54 38];

%% Load OpNav shape

filename='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_150604_GRAVITY_SPC/SHAPE_SPC150604.bds';

step=0.5;
Npts=10000;

[x_rand,y_rand,z_rand]=LoadOpNavShape(filename,Npts,'rand');
[x_grid,y_grid,z_grid]=LoadOpNavShape(filename,step,'grid');


[~,fi_rand,rad_rand]=cart2sph(x_rand,y_rand,z_rand);
fi_rand=fi_rand*180/pi;
cond = ((fi_rand < 75) & (fi_rand > -75));

r_rand=[x_rand(cond)' y_rand(cond)' z_rand(cond)'];



r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);

V = Mesh2Volume(x_grid,y_grid,z_grid);

%% Fit test

[ center, radii, evecs, v ] = ...
    biaxial_ellipsoid_fit(r_rand);

%% Compute the best-fit ellipsoid
[ center, radii, evecs, v ] = ellipsoid_fit( r_rand, 0 );
[ center_xy, radii_xy, evecs_xy, v_xy ] = ellipsoid_fit( r_rand, 2 , 'xy' );
[ center_bi, radii_bi, evecs_bi, v_bi ] = ellipsoid_fit( r_rand, 4);

% [ center, radii, evecs, v ] = ellipsoid_fit_err( r_rand, 4 );
% radii_std=std(radii,0,2);
% center_std=std(center,0,2);

offset=norm(center);
[lambda_ax,fi_ax]=AxesOrientation(evecs);
a = radii(1); b = radii(2); c = radii(3);
ab=sqrt(a*b);
fp = (a - c)/a;
fq = (a - b)/c;
fabp = (ab - c)/ab;

Rvol=(a*b*c)^(1/3);

H_grid = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center,radii,evecs);

H_rand = EllipsoidHeight([x_rand(:) y_rand(:) z_rand(:)],...
    center,radii,evecs);

H_grid_xy = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center_xy,radii_xy,evecs_xy);

H_grid_bi = EllipsoidHeight([x_grid(:) y_grid(:) z_grid(:)],...
    center_bi,radii_bi,evecs_bi);

H_grid=reshape(H_grid,size(x_grid));
H_grid_xy=reshape(H_grid_xy,size(x_grid));
H_grid_bi=reshape(H_grid_bi,size(x_grid));


[lambda,fi,rad]=cart2sph(x_grid,y_grid,z_grid);

%% Write Ellipsoidal elevation grid
WriteXYZ(lambda*180/pi,fi*180/pi,H_grid,'~/Dawn/CeresShapeModel/OpNav/OpNav5/OpNav5_3el.xyz');
WriteXYZ(lambda*180/pi,fi*180/pi,H_grid_xy,'~/Dawn/CeresShapeModel/OpNav/OpNav5/OpNav5_2xy.xyz');
WriteXYZ(lambda*180/pi,fi*180/pi,r_grid,'~/Dawn/CeresShapeModel/OpNav/OpNav5/OpNav5_rad.xyz');
WriteXYZ(lambda*180/pi,fi*180/pi,H_grid_bi,'~/Dawn/CeresShapeModel/OpNav/OpNav5/OpNav5_bi.xyz');

%% Elevation Histogram

figure; hold on;
set(gca,'FontSize',20);

hist(H_rand,20,'b');

xlabel('Elevation [km]','FontSize',20);
ylabel('Number of points []','FontSize',20);

box on;

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

[sdl,l]=plm2spec(lmcosi);
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
set(gcf, 'Units','centimeters', 'Position',[0 0 imsize])
box on;

% xlim([0 10]);

% q=IsotropicRatio(lmcosi(2:end,:),lmcosi(2:end,:));
% 
% figure; hold on;
% plot(1:numel(q),q,'-bo')
% set(gca,'Yscale','log');


%%

[k_Ceres,sdl_Ceres]=PowerSpectrum(lmcosi);

figure; hold on;
set(gca,'FontSize',20);
set(gca,'XScale','log');
set(gca,'YScale','log');
plot(k_Ceres,sdl_Ceres,'-k')











