
im = imread('~/Desktop/hb6.png');
im = fliplr(double(im(:,:,1))')/255;
im = 1-im;

im=im(30:end-30,30:end-30);

x = 1:size(im,1);
y = 1:size(im,2);

[y,x]=meshgrid(y,x);

nx = size(x,1);
ny = size(y,2);

figure; hold on;
pcolor(x,y,im); shading flat;

min_lat=-30;
max_lat=30;

min_lon=0;
max_lon=330;

lon = linspace(min_lon,max_lon,nx);
lat = linspace(min_lat,max_lat,ny);

[lat,lon]=meshgrid(lat,lon);

AGUaxes;
pcolorm(lat,lon,im);

step = 0.05;
lati=min_lat:step:max_lat;
loni=min_lon:step:max_lon;
[lati,loni]=meshgrid(lati,loni);

imi = griddata(lon,lat,im,loni,lati,'nearest');
imi=imi';
 
% imi(isnan(imi))=0;
% 
% AGUaxes;
% pcolorm(lati,loni,imi);
% colorbar
% 
% radius = 2;
% [x,y,z] = sph2cart(loni/180*pi,lati/180*pi,imi+radius);
% 
% figure; hold on;
% surf(x,y,z,imi+radius);
% StandardLight;

% filename='/Users/antonermakov/Dawn/CeresShapeModel/OpNav/OpNav5/ceres_opnav5_512.bds';
% 
% [x_grid,y_grid,z_grid]=LoadOpNavShape(filename,step,'grid');
% [lon_grid,lat_grid,r_grid]=cart2sph(x_grid,y_grid,z_grid);
% 
% r_tot = r_grid+30*fliplr(imi);
% 
% [x_tot,y_tot,z_tot]=sph2cart(lon_grid,lat_grid,r_tot);
% 
% figure; hold on;
% surf(x_tot,y_tot,z_tot,r_tot);
% StandardLight;
% view(180,0)



