ccc

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Initial parameters
G                       = 6.67e-11;
r1                      = 470000;
Rref                    = 470000;
aref                    = 470000;
cref                    = 470000;
GM                      = 62.68e9;
density_contrast        = 0.1;
MaxDegreeTopo           = 200;
MaxTopoPower            = 3;
MaxDegreeGrav           = fix(MaxDegreeTopo/MaxTopoPower);
lmcosi_shape            = plm2rnd(MaxDegreeTopo,-3,1);
lmcosi_shape(1,3)       = abs(lmcosi_shape(1,3));
lmcosi_shape(:,3:4)     = (lmcosi_shape(:,3:4) * r1);
lmcosi_shape(2:end,3:4) = lmcosi_shape(2:end,3:4)/50e2;
lmcosi_shape(2:3,3:4)   = 0;

M = GM/G;

[r,lon,lat] = plm2xyz(lmcosi_shape);
[lon,lat]   = meshgrid(lon,lat);

[x,y,z]     = sph2cart(lon/180*pi,lat/180*pi,r);

%% Compute gravity from topography

[xref,yref,zref]=TriEllRadVec(...
    lat/180*pi,lon/180*pi,aref,aref,cref,'xyz');

lmcosi_gt=Topo2Grav(r,Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

[ax_gt,ay_gt,az_gt]=GravityAcceleration(...
    GM,Rref,lmcosi_gt,xref,yref,zref);

[g_up_gt,~,~]=GravityComponents(...
    ax_gt,ay_gt,az_gt,xref,yref,zref,aref,cref);


%% Isostatic

rho1 = 2000;
rho2 = 2200;

r2 = ((3*M/pi-4*r1^3*rho1)^(1/3))/((4*(rho2-rho1))^(1/3));

D_comp = r1 - r2;

t = FindCrustalRoot(r1,D_comp,r-r1,rho1,rho2-rho1);

r2_grid = TriEllRadVec(lat/180*pi,lon/180*pi,r1-D_comp,r1-D_comp,r1-D_comp,'rad');
r2_grid = r2_grid - t;

% r2_grid = r_grid - D_comp - h;

[x2,y2,z2] = sph2cart(lon/180*pi,lat/180*pi,r2_grid);

lmcosi_gt2_isos=Topo2Grav(r2_grid,Rref,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

V1 = Mesh2Volume(x,y,z);
V2 = Mesh2Volume(x2,y2,z2);

rho_mean = M/V1;

M1_isos = rho1*V1;
M2_isos = (rho2-rho1)*V2;

w_isos = [M1_isos/(M1_isos+M2_isos) M2_isos/(M1_isos+M2_isos)];

lmcosi_gt_isos = WeightSumExpansion(w_isos,{lmcosi_gt,lmcosi_gt2_isos});
%% Compute isostatic gravity

[ax_gt_isos,ay_gt_isos,az_gt_isos]=GravityAcceleration(...
    GM,Rref,lmcosi_gt_isos,xref,yref,zref);

[g_up_gt_isos,~,~]=GravityComponents(...
    ax_gt_isos,ay_gt_isos,az_gt_isos,xref,yref,zref,aref,cref);

%% Plot maps

figure; hold on;
pcolor(lon,lat,r); shading flat; colorbar
title('shape');

figure; hold on;
pcolor(lon,lat,g_up_gt); shading flat; colorbar
title('Homo gravity');

figure; hold on;
pcolor(lon,lat,g_up_gt_isos); shading flat; colorbar
title('Isostatic gravity');

%% Compute Admittance

[n,Z_gt] = SphericalHarmonicAdmittance(lmcosi_gt,lmcosi_shape,GM,Rref);
[n,Z_isos] = SphericalHarmonicAdmittance(lmcosi_gt_isos,lmcosi_shape,GM,Rref);

R0    = lmcosi_shape(1,3);
Z_theo  = 3./(2*(n)+1).*((R0/Rref).^n).*(n+1)*GM/(Rref^2)*1e5/R0*1000;

Z_theo_isos  = 3./(2*(n)+1).*rho1/rho_mean.*(1-(1-D_comp./r1).^n).*...
    ((R0/Rref).^n).*(n+1)*GM/(Rref^2)*1e5/R0*1000;

figure;
hold on;box on;grid on;
set(gcf, 'Units','centimeters', 'Position',im_size);
set(gcf, 'PaperPositionMode','auto');
set(gca, 'FontSize',fntsize);

% h_gt    = plot(n,Z_gt,'ro-','LineWidth',3);
h_theo = plot(n,Z_theo,'ro--','LineWidth',3);
% plot(n,Z_isos,'bo-','LineWidth',3);
% plot(n,Z_theo_isos,'bo--','LineWidth',3);

xlabel('Degree','FontSize',12);
ylabel('Admittance [mGal/km]','FontSize',12);



%%
load '2layer_solution.mat';

rho1_Jh_grid = 1000:100:2000;
cond = (rho1_Jh_grid < min(rho1_Jh)) |  (rho1_Jh_grid > max(rho1_Jh));
rho1_Jh_grid(cond) = [];

rho2_Jh_grid = interp1(rho1_Jh(~isnan(rho1_Jh)),rho2_Jh(~isnan(rho1_Jh)),rho1_Jh_grid);
r2_Jh_grid   = interp1(rho1_Jh(~isnan(rho1_Jh)),r2_Jh(~isnan(rho1_Jh)),rho1_Jh_grid);

for i = 1:numel(rho1_Jh_grid)
    D_comp = r1 - r2_Jh_grid(i);
    Z_theo_isos  = 3./(2*(n)+1).*rho1_Jh_grid(i)/rho_mean.*(1-(1-D_comp./r1).^n).*...
    ((R0/Rref).^n).*(n+1)*GM/(Rref^2)*1e5/R0*1000;
   
    h_isos = plot(n,Z_theo_isos,'m--');
    text(n(1),Z_theo_isos(1),num2str(rho1_Jh_grid(i)));
    
end

load Z_g_obs.mat;
h_obs = plot(n,Z_g,'-ok','LineWidth',3);

xlim([2 5]);
ylim([-60 120]);

set(gca,'XTick',2:5);

legend([h_obs h_isos h_theo],{'Observed','Isostatic linear','Homo','Lin homo'},'FontSize',12);







