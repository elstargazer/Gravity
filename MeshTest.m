ccc

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';
[~,shapename,~] = fileparts(shape_filename) ;
shape_full_filename = [shape_folder shape_filename];

filename_mesh = '../CeresFE/FE/mesh_gibbon.inp';

%% Read Ceres shape 
step = 1;
[x_grid,y_grid,z_grid]=ReadSPC(shape_full_filename,step,'grid');

x_grid = x_grid*1000;
y_grid = y_grid*1000;
z_grid = z_grid*1000;

r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);

figure; hold on;
surf(x_grid,y_grid,z_grid,r_grid); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z'); 

lmcosi_real_shape = xyz2plm(flipud(r_grid'));

[sdl,l] = plm2spec(lmcosi_real_shape);
p = polyfit(log10(l(21:end)),log10(sdl(21:end)),1);

figure; hold on;
set(gca,'XScale','log');
set(gca,'YScale','log');

plot(l,sdl);

xlabel('Degree','FontSize',20);
ylabel('Power','FontSize',20);

%% Generate Random mesh

r_sphe  = 470000;
r_core  = 300000;
r_cube  = 100000;
n_cube  = 10;
n_core  = 10;
n_sphe  = 10;

L          = 150;
beta       = p(1);
intercept  = 10^p(2);

% Shape2Mesh_gibbon(r_sphe, r_core, n_core, n_sphe, ...
%     shape_full_filename, filename_mesh);

RandomShape2Mesh_gibbon(r_sphe, r_cube, n_cube, n_sphe, ...
    filename_mesh,L, beta, intercept)

%% Generate two-layer random mesh
% 
% % generate outer surface
% L    = 10;
% beta = -2;
% 
% lmcosi1=plm2rnd(L,beta,1,3);
% 
% lmcosi1(1,3) = abs(lmcosi1(1,3)*r_sphe);
% lmcosi1(2:end,3:4) = lmcosi1(2:end,3:4)*2e4;
% 
% % generate inner surface 
% L    = 10;
% beta = -3;
% 
% lmcosi2=plm2rnd(L,beta,1,3);
% 
% lmcosi2(1,3) = abs(lmcosi2(1,3)*r_core);
% lmcosi2(2:end,3:4) = lmcosi2(2:end,3:4)*0.5e4;
% 
% 
% figure; hold on;
% axis equal;
% alpha(0.5);
% 
% [r1,lon,lat] = plm2xyz(lmcosi1);
% [lon,lat] = meshgrid(lon,lat);
% [x1,y1,z1] = sph2cart(lon/180*pi,lat/180*pi,r1);
% 
% surf(x1,y1,z1);
% 
% 
% [r2,lon,lat] = plm2xyz(lmcosi2);
% [lon,lat]=meshgrid(lon,lat);
% [x2,y2,z2] = sph2cart(lon/180*pi,lat/180*pi,r2);
% 
% surf(x2,y2,z2);
% 
% 
% alpha(0.5);
% 
% TwoLayer2Mesh_gibbon(lmcosi1, lmcosi2, n_sphe, n_core, n_cube, r_cube, filename_mesh);

%% Read outer surface and compute power specturm

% L = 60;
% filename = '../CeresFE/FE/outer_points.txt';
% [sdl,l] = file2spec(filename,L);

% figure; hold on;
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% 
% plot(l,sdl,'r-');
% 
% xlabel('Degree','FontSize',20);
% ylabel('Power','FontSize',20);

%% Generate Jacobi ellipsoid mesh
% T   = 4.6;
% rho = 2169;
% [fh,fval] = HydrostaticStateExact3Ax(r_sphe,T,rho,0.05,0.025);
% [a,b,c] = fr2abc(r_sphe,fh(1),fh(2));
% Ellipsoid2Mesh_gibbon([a b c], r_core, n_core, n_sphe, filename_mesh)

% nrow = 10;
% Shape2Mesh_cubesphere(shape_full_filename,filename_mesh,nrow);
 