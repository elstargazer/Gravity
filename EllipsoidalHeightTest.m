%% Load shape model


ShapeModelFileName='~/Dawn/Balmino/VestaTest/VestaPOSTVESTAshape_geoc_elev12m.grd';

[lambdai,fii,ri]=ReadGRD(ShapeModelFileName);

ri=ri*1000;
lambdai=lambdai/180*pi;
fii=fii/180*pi;

[xi,yi,zi]=sph2cart(lambdai,fii,ri);

figure
surf(xi,yi,zi,ri);
StandardLight

%% Ellipsoid parameters

a=280000;
b=





[min_dist, f_min] = distanceEllipsePoints(XYZ, a,b,c,u,v)