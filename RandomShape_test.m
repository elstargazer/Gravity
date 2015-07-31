% ccc;

%% Initial parameters

r_sph = 470000;
cell_type = 'quad';

init_mesh_filename = '../CeresFE/FE/mesh_test.inp';
[path,name,ext] = fileparts(init_mesh_filename);
deformed_mesh_filename = [path '/' name '_def' ext];
deformed_mesh_quad_filename = [path '/' name '_def_quad' ext];

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';

filename_quad = '../CeresFE/FE/outer_points.txt';


%% Read init sphere mesh

meshStruct = Read_ucd(init_mesh_filename);

x = meshStruct.V(:,1);
z = meshStruct.V(:,2);

figure; hold on;
plot(x,z,'.');
axis equal;

real_shape_full_filename = [shape_folder shape_filename];

step = 1;
Npts = 100000;
[x_grid,y_grid,z_grid]=ReadSPC(real_shape_full_filename,step,'grid');

x_grid = x_grid*1000;
y_grid = y_grid*1000;
z_grid = z_grid*1000;

r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);

figure; hold on;
surf(x_grid,y_grid,z_grid,r_grid); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z'); 

lmcosi_real_shape = xyz2plm(flipud(r_grid'));
[r_sh, lon_sh, lat_sh] = plm2xyz(lmcosi_real_shape);
[lon_sh,lat_sh] = meshgrid(lon_sh,lat_sh);
[x_sh,y_sh,z_sh] = sph2cart(lon_sh/180*pi,lat_sh/180*pi,r_sh);

figure; hold on;
surf(x_sh,y_sh,z_sh,r_sh); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z'); 

%% Ceres spectrum

[sdl,l] = plm2spec(lmcosi_real_shape);

plot_spec = figure; hold on;
set(gca,'FontSize',20);
set(gca,'XScale','log');
set(gca,'YScale','log');
plot(l,sdl);
box on;
xlabel('Degree','FontSize',20);
ylabel('Power [km^{2}]','FontSize',20);

gi = ginput(2);

minn=round(gi(1,1));
maxn=round(gi(2,1));

cond = ((l>minn) & (l<maxn));

l_c   = l(cond);
sdl_c = sdl(cond);

p = polyfit(log10(l_c),log10(sdl_c),1);

ipt1 = polyval(p,log10(1));

%% Generate random spectrum

L    = 180;
% beta = p(1);
% 
% [lmcosi_shape,bta,bto,sdl,el]=plm2rnd(L,beta,1);
% lmcosi_shape(2:end,3:4) = lmcosi_shape(2:end,3:4)*sqrt((10^ipt1));
% lmcosi_shape(1,3) = lmcosi_real_shape(1,3);
% 
% lmcosi_shape(2,3:4) = 0;
% lmcosi_shape(3,3:4) = 0;
% 
% [sdl_r,l_r] = plm2spec(lmcosi_shape);

% lmcosi_shape(2:end,3:4) = 0;
lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = lmcosi_real_shape(1,3);
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;


for n=2:2:lmcosi_shape(end,1)  
%     lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*sdl_r(n+1));
    lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*10^polyval(p,log10(n)));
end

[sdl_r,l_r] = plm2spec(lmcosi_shape);

figure(plot_spec);
plot(l_r,sdl_r,'.g','MarkerSize',15);


%% Compute limb

[lon,lat,r] = cart2sph(x,0,z);
r_shape = plm2xyz(lmcosi_shape,lat*180/pi,lon*180/pi);

r_new = r.*r_shape;
[x_new, y_new, z_new] = sph2cart(lon,lat,r_new);

figure; hold on;
plot(x_new,z_new,'.');


meshStruct_def.E = meshStruct.E;
meshStruct_def.V = [x_new z_new zeros(size(x_new))];

xv1 = meshStruct_def.V(meshStruct_def.E(:,1),1);
xv2 = meshStruct_def.V(meshStruct_def.E(:,2),1);
xv3 = meshStruct_def.V(meshStruct_def.E(:,3),1);
xv4 = meshStruct_def.V(meshStruct_def.E(:,4),1);

zv1 = meshStruct_def.V(meshStruct_def.E(:,1),2);
zv2 = meshStruct_def.V(meshStruct_def.E(:,2),2);
zv3 = meshStruct_def.V(meshStruct_def.E(:,3),2);
zv4 = meshStruct_def.V(meshStruct_def.E(:,4),2);

eps = 1e-5;
ind_quad = find((xv1 > -eps) & (xv2 > -eps) & (xv3 > -eps) & (xv4 > -eps) & ...
    (zv1 > -eps) & (zv2 > -eps) & (zv3 > -eps) & (zv4 > -eps));
    
meshStruct_def_quad.E = meshStruct_def.E(ind_quad,:);
meshStruct_def_quad.V = meshStruct_def.V;

Write_ucd(meshStruct_def_quad,deformed_mesh_quad_filename,cell_type)
Write_ucd(meshStruct_def,deformed_mesh_filename,cell_type)

axis equal;

%% limb spectrum
L = 40;

lmcosi_limb = quad2plm(filename_quad,L);

[sdl_limb,l_limb] = plm2spec(lmcosi_limb);
figure(plot_spec);
plot(l_limb,sdl_limb,'.m','MarkerSize',15);

[r_limb_sh, lon_limb_sh, lat_limb_sh] = ...
    plm2xyz(lmcosi_limb);

[lon_limb_sh,lat_limb_sh] = meshgrid(lon_limb_sh,lat_limb_sh);
[x_limb_sh,y_limb_sh,z_limb_sh] = ...
    sph2cart(lon_limb_sh/180*pi,lat_limb_sh/180*pi,r_limb_sh);

figure; hold on;
surf(x_limb_sh,y_limb_sh,z_limb_sh,r_limb_sh); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z'); 






