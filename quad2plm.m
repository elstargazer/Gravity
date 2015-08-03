function lmcosi_limb = quad2plm(filename,L)

data = load(filename);

x_limb = data(:,1);
z_limb = data(:,2);

eps = 1;
cond_bdry = (abs(x_limb) < eps | abs(z_limb < eps));

x_limb(cond_bdry) = [];
z_limb(cond_bdry) = [];
 
% figure; hold on;
% plot(x_limb,z_limb,'.');

x_limb = [x_limb; -x_limb; -x_limb;  x_limb];
z_limb = [z_limb;  z_limb; -z_limb; -z_limb];

[~,lat,r_limb] = cart2sph(x_limb,0,z_limb);

% lat_interp = (-90:1:90)/180*pi;
% r_limb_interp = interp1(lat,r_limb,lat_interp);

lon = linspace(-pi,pi,2*numel(x_limb));
[lon,lat] = meshgrid(lon,lat);
r_limb = repmat(r_limb,[1 2*numel(x_limb)]);

lmcosi_limb = xyz2plm(r_limb(:),L,'irr',lat(:)*180/pi,lon(:)*180/pi);

% [r_limb_sh, lon_limb_sh, lat_limb_sh] = ...
%     plm2xyz(lmcosi_limb);
% 
% [lon_limb_sh,lat_limb_sh] = meshgrid(lon_limb_sh,lat_limb_sh);
% [x_limb_sh,y_limb_sh,z_limb_sh] = ...
%     sph2cart(lon_limb_sh/180*pi,lat_limb_sh/180*pi,r_limb_sh);
% 
% figure; hold on;
% surf(x_limb_sh,y_limb_sh,z_limb_sh,r_limb_sh); shading interp
% axis equal
% xlabel('x'); ylabel('y'); zlabel('z'); 

% figure;
% scatter(lon(:),lat(:),10,r_limb(:),'filled');
% 
% 
% [x_limb,y_limb,z_limb] = sph2cart(lon,lat,r_limb);
% figure; hold on;
% axis equal;
% plot3(x_limb(:),y_limb(:),z_limb(:),'.'); shading interp;
% StandardLight;
% 
% % ellipsoid
% ang = linspace(0,2*pi,100);
% ell = [479649.298911 449096.692424];
% xfit_ell = ell(1)*cos(ang);
% zfit_ell = ell(2)*sin(ang);
% figure; hold on;
% plot(x_limb,z_limb,'r.');
% plot(xfit_ell,zfit_ell,'b.');



