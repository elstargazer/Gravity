% ccc

aref = 500000;
cref = 500000;
rho = 2161;

%% Load gravity models

NTess = 4;
shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150828_512.bds';
fig_folder='~/Dawn/Papers/CeresPaper1/';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];
[x1,y1,z1,FV1] = DSK2TRI(full_filename,NTess);
x1=x1*1000;
y1=y1*1000;
z1=z1*1000;

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPG/Survey/';
shape_filename='global.bds';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];
[x2,y2,z2,FV2] = DSK2TRI(full_filename,NTess);

step = 1;
lat=(-90:step:90) * cspice_rpd();
lon=(0:step:360)  * cspice_rpd();
[lati, loni] = meshgrid(lat, lon);

%% Find min ellipsoid reference surface
[A , c] = MinVolEllipse([x2(:)'; y2(:)'; z2(:)'], 1e-4);
[U, Q, V] = svd(A);
ell_min=1./sqrt(diag(Q));

aref = ell_min(3);
bref = ell_min(2);
cref = ell_min(1);

[xref,yref,zref]=TriEllRadVec(lati,loni,aref,aref,cref,'xyz');

%% Compute gravity 
s=size(FV1);
rho=ones(1,s(1))*rho;

[ax1,ay1,az1]=GravityAccelerationTriDen(x1,y1,z1,FV1,xref,yref,zref,rho);
[ax1_up,ax1_east,ax1_north]=GravityComponents(ax1,ay1,az1,xref,yref,zref,aref,cref);

[ax2,ay2,az2]=GravityAccelerationTriDen(x2,y2,z2,FV2,xref,yref,zref,rho);
[ax2_up,ax2_east,ax2_north]=GravityComponents(ax2,ay2,az2,xref,yref,zref,aref,cref);

 dg = sqrt((ax1-ax2).^2+(ay1-ay2).^2+(az1-az2).^2)*1e5;
 dg_up = ax1_up - ax2_up;

AGUaxes;
pcolorm(lati*180/pi,loni*180/pi,dg_up*1e5); 
cbar = colorbar('FontSize',20);
ylabel(cbar,'Acceleration difference [mGal]','FontSize',20);






