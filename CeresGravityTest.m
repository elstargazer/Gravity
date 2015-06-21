ccc

filename='/Users/antonermakov/Dawn/CeresShapeModel/OpNav/OpNav5/ceres_opnav5_512.bds';
filename_grd='/Users/antonermakov/Dawn/CeresShapeModel/OpNav/OpNav5/ceres_opnav5_512.bds';

step=0.5;
T=9.074;
omega=2*pi/T/3600;

mu=62.68e9;
rho=2161;
Rref_grav=500000;
MaxDegreeTopo=360;
MaxDegreeGrav=9;
MaxTopoPower=5;

[x_grid,y_grid,z_grid]=LoadOpNavShape(filename,step,'grid');

x_grid=x_grid*1000;
y_grid=y_grid*1000;
z_grid=z_grid*1000;

[lambda_grid, fi_grid, r_grid]=cart2sph(x_grid,y_grid,z_grid);

lmcosi_grav=TopoSH2GravitySH(r_grid',mu,rho,Rref_grav,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

ag=520048.2095350715;
cg=487664.3599076933;


[xs,ys,zs]=TriEllRadVec(fi_grid',lambda_grid',ag,ag,cg,'xyz');

[ax,ay,az]=GravityAcceleration(mu,Rref_grav,lmcosi_grav,xs,ys,zs);
[g_up,g_east,g_north]=GravityComponents(ax,ay,az,xs,ys,zs,ag,cg);

lmcosi_ell=SHRotationalEllipsoid(ag,cg,10,Rref_grav);

[ax_norm,ay_norm,az_norm]=GravityAcceleration(mu,Rref_grav,lmcosi_ell,xs,ys,zs);
[g_up_norm,g_east_norm,g_north_norm]=GravityComponents(ax_norm,ay_norm,az_norm,...
    xs,ys,zs,ag,cg);
% 
% AGUaxes
% pcolorm(fi_grid'*180/pi,lambda_grid'*180/pi,1e5*(g_up));
% cbar=colorbar('FontSize',25);
% ylabel(cbar,'Gravity Acceleration','FontSize',25);
% 
% AGUaxes
% pcolorm(fi_grid'*180/pi,lambda_grid'*180/pi,1e5*(g_up_norm));
% cbar=colorbar('FontSize',25);
% ylabel(cbar,'Gravity Acceleration','FontSize',25);

fa=(1e5*(g_up-g_up_norm));

AGUaxes
pcolorm(fi_grid'*180/pi,(lambda_grid'*180/pi),fa);
cbar=colorbar('FontSize',25);
ylabel(cbar,'Free-Air anomaly [mGal]','FontSize',25);


a=482788.0;
c=447800.0;
[B,L,H]=XYZ2BLH(x_grid,y_grid,z_grid,a,Eccentricity(a,c));

minH=min(H(:));
maxH=max(H(:)); 
Hm=(H-minH)/(maxH-minH);

AGUaxes
surfm(fi_grid'*180/pi,lambda_grid'*180/pi,H'/1000,Hm'-2);
cbar=colorbar('FontSize',25);
ylabel(cbar,'Ellipsoidal Height [km]','FontSize',25);
zlim([0 1]-2)

AGUaxes
surfm(fi_grid'*180/pi,lambda_grid'*180/pi,fa,Hm'-2);
cbar=colorbar('FontSize',25);
ylabel(cbar,'Free-Air anomaly [mGal]','FontSize',25);
zlim([0 1]-2)

lighting phong
shading interp
light('Position',[-1 1 1],'Style','infinite');


% AGUaxes
% pcolorm(fi_grid'*180/pi,lambda_grid'*180/pi,H_grid_bi');
% cbar=colorbar('FontSize',25);
% ylabel(cbar,'Gravity Acceleration','FontSize',25);














