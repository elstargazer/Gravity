function H = EllipsoidHeight(rvec,center,radii,evecs)

x_raw=rvec(:,1);
y_raw=rvec(:,2);
z_raw=rvec(:,3);

x_mod=x_raw-center(1);
y_mod=y_raw-center(2);
z_mod=z_raw-center(3);

r_mod=[x_mod y_mod z_mod]*evecs;

x_mod=r_mod(:,1);
y_mod=r_mod(:,2);
z_mod=r_mod(:,3);

[ center_mod, radii_mod, evecs_mod, v_mod ] = ellipsoid_fit( [x_mod y_mod z_mod], 0 );

[center_mod center]
[radii_mod radii]
[evecs_mod evecs]

[ ~, H ] = cspice_nearpt( r_mod', radii(1), radii(2), radii(3) );
