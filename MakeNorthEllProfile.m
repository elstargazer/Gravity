

MaxDegreeTopo=500;
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);



Resolution=0.2;
 
[ri_shape,lambda,fi,~]=plm2xyz(lmcosi_shape,Resolution);

[lambdai,fii]=meshgrid(lambda/180*pi,fi/180*pi);

close all
%% Lambda = 0

lambdap=0/180*pi;

Cond=((lambdai==lambdap) | (lambdai==lambdap+pi));

[x,y,z]=sph2cart(lambdai(Cond),fii(Cond),ri_shape(Cond));

[~,~,H]=XYZ2BLH(x,y,z,280917.27,Eccentricity(280917.27,226248.16));

vesta_fig=figure; hold on;

p = patch( isosurface( x_fit, y_fit*0, z_fit, Ellipsoid, 1 ) );
set( p, 'FaceColor', 'b', 'EdgeColor', 'none' );

alpha(0.4)

set(gca,'FontSize',20);
plot3(x/1000,y/1000,z/1000,'k','LineWidth',2);

view(lambdap*180/pi,0)

xlabel('x [km]','FontSize',20);
ylabel('y [km]','FontSize',20);
zlabel('z [km]','FontSize',20);

box on;
axis equal
view(lambdap*180/pi,0)

xlim([-300 300])
ylim([-300 300])
zlim([-300 300])

set(gcf, 'Units','centimeters', 'Position',[0 0 17 17])
set(gcf, 'PaperPositionMode','auto')

print(vesta_fig, '-dpsc2', 'VestaProfileNorthEllLambda0.eps');

%% Lambda = 90

lambdap=90/180*pi;
Cond=((lambdai==lambdap) | (lambdai==lambdap+pi));
[x,y,z]=sph2cart(lambdai(Cond),fii(Cond),ri_shape(Cond));

vesta_fig=figure; hold on;
p = patch( isosurface( x_fit*0, y_fit, z_fit, Ellipsoid, 1 ) );
set( p, 'FaceColor', 'b', 'EdgeColor', 'none' );

alpha(0.4)
set(gca,'FontSize',20);
plot3(x/1000,y/1000,z/1000,'k','LineWidth',2);

view(lambdap*180/pi,0)

xlim([-300 300])
ylim([-300 300])
zlim([-300 300])

xlabel('x [km]','FontSize',20);
ylabel('y [km]','FontSize',20);
zlabel('z [km]','FontSize',20);

box on;
axis equal
view(lambdap*180/pi,0)

xlim([-300 300])
ylim([-300 300])
zlim([-300 300])

set(gcf, 'Units','centimeters', 'Position',[0 0 17 17])
set(gcf, 'PaperPositionMode','auto')

print(vesta_fig, '-dpsc2', 'VestaProfileNorthEllLambda90.eps');