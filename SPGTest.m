ccc

G=6.67e-11;
GM=62.8981e9;

filename='~/Dawn/CeresShapeModel/SPG/PreSurvey/global.xyz';
data=load(filename);

cubfile='~/Dawn/CeresShapeModel/SPG/PreSurvey/global.cub';
[label, core, isis_version] = read_isis_cub(cubfile);

core=core(2:end-1,:);
core(:,1)=[];

core=flipud(core);

long=linspace(0,360,size(core,2))/180*pi;
latg=linspace(-90,90,size(core,1))/180*pi;

[long,latg]=meshgrid(long,latg);

[xg,yg,zg]=sph2cart(long,latg,core);

sg = size(xg);

xg=xg(:);
yg=yg(:);
zg=zg(:);


% WriteXYZ(data(:,2),data(:,1),data(:,3),'../CeresShapeModel/SPC/CeresShapeFeb16.xyz');
% WriteXYZ(lon*180/pi,lat*180/pi,rad,'../CeresShapeModel/SPG/2015_02_26.Ceres-RC2.SPG.2ppd.gaps.sph');

% lat=data(:,1)/180*pi;
% lon=data(:,2)/180*pi;
% rad=data(:,3);
% 
% [x_raw,y_raw,z_raw]=sph2cart(lon,lat,rad);

x_raw=data(:,1);
y_raw=data(:,2);
z_raw=data(:,3);

r_raw=[x_raw y_raw z_raw];
rg = [xg(:) yg(:) zg(:)];

[lon,lat,rad]=cart2sph(data(:,1),data(:,2),data(:,3));

%% SH analysis
L=360;

% lmcosi_topo=xyz2plm(rad,L,'irr',lat*180/pi,lon*180/pi);
lmcosi_topo=xyz2plm(core,L);

[r_sh,lon_sh,lat_sh]=plm2xyz(lmcosi_topo,1);

[lon_sh,lat_sh]=meshgrid(lon_sh,lat_sh);
lon_sh=lon_sh/180*pi;
lat_sh=lat_sh/180*pi;

[x_sh,y_sh,z_sh]=sph2cart(lon_sh,lat_sh,r_sh);

figure; hold on;
surf(x_sh,y_sh,z_sh,r_sh); shading interp;
StandardLight;

% figure;
% plot3(x_raw,y_raw,z_raw,'.');

tic
[ center, radii, evecs, v ] = ellipsoid_fit( [x_raw y_raw z_raw], 0 );
toc

[ center4, radii4, evecs4, v4 ] = ellipsoid_fit( [x_raw y_raw z_raw], 4 );
[ center5, radii5, evecs5, v5 ] = ellipsoid_fit( [x_raw y_raw z_raw], 5 );

%WriteXYZ(lon*180/pi,lat*180/pi,rad,'~/Dawn/CeresShapeModel/SPG/PreSurvey/global.');

%% Modifying to principal ellipsoid
x_mod=x_raw-center(1);
y_mod=y_raw-center(2);
z_mod=z_raw-center(3);

xg_mod = xg-center(1);
yg_mod = yg-center(2);
zg_mod = xg-center(3);

r_mod=[x_mod y_mod z_mod]*evecs;
rg_mod=[xg_mod yg_mod zg_mod]*evecs;

x_mod=r_mod(:,1);
y_mod=r_mod(:,2);
z_mod=r_mod(:,3);

xg_mod=rg_mod(:,1);
yg_mod=rg_mod(:,2);
zg_mod=rg_mod(:,3);

[ center_mod, radii_mod, evecs_mod, v_mod ] = ellipsoid_fit( [x_mod y_mod z_mod], 0 );

[ ~, H ] = cspice_nearpt( r_mod', radii_mod(1), radii_mod(2), radii_mod(3) );
[ ~, Hg ] = cspice_nearpt( rg_mod', radii_mod(1), radii_mod(2), radii_mod(3) );

WriteXYZ(long(:)*180/pi,latg(:)*180/pi,Hg/1000,'~/Dawn/CeresShapeModel/SPG/PreSurvey/global_H.dat');


%% Fititng biaxial ellipsoid

[ center_bi, radii_bi, evecs_bi, v_bi] = ellipsoid_fit( [x_raw y_raw z_raw], 2, 'xy' );

x_bi=x_raw-center_bi(1);
y_bi=y_raw-center_bi(2);
z_bi=z_raw-center_bi(3);

xg_bi=xg-center_bi(1);
yg_bi=yg-center_bi(2);
zg_bi=zg-center_bi(3);

r_bi=[x_bi y_bi z_bi];
rg_bi=[xg_bi yg_bi zg_bi];

a_bi=radii_bi(1);
c_bi=radii_bi(3);

f_bi=(a_bi-c_bi)/a_bi

[ ~, H_bi ] = cspice_nearpt( r_bi', ...
    radii_bi(1), radii_bi(2), radii_bi(3) );

[ ~, Hg_bi ] = cspice_nearpt( rg_bi', ...
    radii_bi(1), radii_bi(2), radii_bi(3) );

WriteXYZ(long(:)*180/pi,latg(:)*180/pi,Hg_bi/1000,...
    '~/Dawn/CeresShapeModel/SPG/PreSurvey/global_Hbi.dat');

%% 

[ ~, H_raw ] = cspice_nearpt( r_raw', 481000, 481000, 447000 );
%WriteXYZ(lon*180/pi,lat*180/pi,H_raw,'../CeresShapeModel/SPG/2015_02_26.Ceres-RC2.SPG.2ppd.gaps.raw');

%% Histograms

figure; hold on;

hist(H_bi,20,'r');
hist(H,20,'b');


%% Ellipsoid parameters

radii_rot=[sqrt(radii(1)*radii(2)) sqrt(radii(1)*radii(2)) radii(3)];

f=(radii_rot(1)-radii_rot(3))/radii_rot(1);
 
a_vec=evecs(:,1);
b_vec=evecs(:,2);
c_vec=evecs(:,3);

a=radii(1);
b=radii(2);
c=radii(3);

ab=sqrt(a*b);

fab_obs=(ab-c)/ab;

fp_obs=(a-c)/a;
fq_obs=(a-b)/a;

rvol=(a*b*c)^(1/3);
vol=4/3*pi*rvol^3;
M=GM/G;
rho_mean=M/vol

[a b c]'/1000
[fab_obs fp_obs fq_obs]'
center/1000

%%

%% Draw fit
% maxd = max( radii );
% step = maxd / 50;
% 
% 
% [ x_fit, y_fit, z_fit ] = meshgrid( -maxd:step:maxd + center(1), ...
%     -maxd:step:maxd + center(2), -maxd:step:maxd + center(3) );
% 
% Ellipsoid = v(1) *x_fit.*x_fit +   v(2) * y_fit.*y_fit + v(3) * z_fit.*z_fit + ...
%           2*v(4) *x_fit.*y_fit + 2*v(5)*x_fit.*z_fit + 2*v(6) * y_fit.*z_fit + ...
%           2*v(7) *x_fit    + 2*v(8)*y_fit    + 2*v(9) * z_fit;
%       
% p = patch( isosurface( x_fit, y_fit, z_fit, Ellipsoid, 1 ) );
% set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
% view( -70, 40 );
% axis vis3d;
% 
% 
% view(0,0)
% 
% zlim([-600000 600000]);
% xlim([-600000 600000]);
% ylim([-600000 600000]);
%  

% plotplm(lmcosi_shape,[],[],2,1,[],[],[]);    

% hold on

% vec_lenght=700000;
 
%% Ellipsoid axis for the North shape
% line([center(1) center(1)+vec_lenght*a_vec(1)],...
%     [center(2) center(2)+vec_lenght*a_vec(2)],...
%     [center(3) center(3)+vec_lenght*a_vec(3)],'Color','r','LineWidth',3,'LineStyle','-');

% line([center_cut(1) center_cut(1)+vec_lenght*b_vec_cut(1)],...
%     [center_cut(2) center_cut(2)+vec_lenght*b_vec_cut(2)],...
%     [center_cut(3) center_cut(3)+vec_lenght*b_vec_cut(3)],'Color','g','LineWidth',3,'LineStyle','-');
% 
% line([center_cut(1) center_cut(1)+vec_lenght*c_vec_cut(1)],...
%     [center_cut(2) center_cut(2)+vec_lenght*c_vec_cut(2)],...
%     [center_cut(3) center_cut(3)+vec_lenght*c_vec_cut(3)],'Color','b','LineWidth',3,'LineStyle','-');

%% Ellipsoid axis for the Total shape

% line([center(1) center(1)+vec_lenght*a_vec(1)],...
%     [center(2) center(2)+vec_lenght*a_vec(2)],...
%     [center(3) center(3)+vec_lenght*a_vec(3)],'Color','r','LineWidth',3,'LineStyle','--');

% line([center(1) center(1)+vec_lenght*b_vec(1)],...
%     [center(2) center(2)+vec_lenght*b_vec(2)],...
%     [center(3) center(3)+vec_lenght*b_vec(3)],'Color','g','LineWidth',3,'LineStyle','--');
% 
% line([center(1) center(1)+vec_lenght*c_vec(1)],...
%     [center(2) center(2)+vec_lenght*c_vec(2)],...
%     [center(3) center(3)+vec_lenght*c_vec(3)],'Color','b','LineWidth',3,'LineStyle','--');

%% Rotation axis
% 
% line([0 0],...
%     [0 0],...
%     [-vec_lenght vec_lenght],'Color','k','LineWidth',4,'LineStyle','-');

% line([0 0],...
%     [-vec_lenght vec_lenght],...
%     [0 0],'Color','k','LineWidth',4,'LineStyle','-');
% 
% line([-vec_lenght vec_lenght],...
%     [0 0],...
%     [0 0],'Color','k','LineWidth',4,'LineStyle','-');

% lambda_eq=0:.1:360;
% fi_eq=lambda_eq.*0;
% 
% 
% [ri_eq,lambda_eq,fi_eq,~]=plm2xyz(lmcosi_shape,fi_eq,lambda_eq);
% 
% [x_eq,y_eq,z_eq]=sph2cart(lambda_eq/180*pi,fi_eq/180*pi,ri_eq');
% 
%  plot3(x_eq/1000,y_eq/1000,z_eq/1000,'-k','LineWidth',4);

% alpha(0.7)
% 
% lambda_a=atan2(a_vec(2),a_vec(1))*180/pi
% lambda_b=atan2(b_vec(2),b_vec(1))*180/pi
% lambda_c=atan2(c_vec(2),c_vec(1))*180/pi
% 
% fi_a=atan2(a_vec(3),norm([a_vec(1) a_vec(2)]))*180/pi
% fi_b=atan2(b_vec(3),norm([b_vec(1) b_vec(2)]))*180/pi
% fi_c=atan2(c_vec(3),norm([c_vec(1) c_vec(2)]))*180/pi

%% 
% 
% [lambdai1,fii1,ri1]=ReadGRD('../CeresShapeModel/SPC/CeresShapeFeb16-2_1d.grd');
% 
% [x1,y1,z1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);
% 
% [~,~,H1]=XYZ2BLH(x1-center(1),y1-center(2),z1-center(3),ab,Eccentricity(ab,c));
% 
% % surf(x1,y1,z1,H1);
% 
% WriteXYZ(lambdai1,fii1,H1,'../CeresShapeModel/SPC/CeresShapeFeb16-2_Ell.xyz');
% 
% 
% figure; hold on;
% pcolor(lambdai1,fii1,H1); shading interp
% caxis([-15000 15000]);
% colorbar


%% Elevation histogram
% 
% figure; hold on;
% set(gca,'FontSize',20);
% 
% hist(H/1000,40)
% xlabel('Elevation [km]','FontSize',20);
% ylabel('Number of points []','FontSize',20);
% box on;




