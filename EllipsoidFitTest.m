[fii,lambdai]=GenerateRandomSphericalCoord(100000);

a0=446;
b0=448;
c0=430;

[x,y,z] = TriEllRadVec(fii,lambdai,a0,b0,c0,'xyz');

plot3(x,y,z,'.');

[ center, radii, evecs, v ] = ...
    triaxial_ellipsoid_fit([x y z]);

[ center2, radii2, evecs2, v ] = ...
    ellipsoid_fit([x y z], 0 );

radii