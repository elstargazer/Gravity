function alt_back=BestFitEllipsoidHeight(x,y,z,center_ell,radii_ell,evecs_ell)


s=size(x);

point_new=[x(:), y(:), z(:)]';


%% Transformation back


point_back(1,:)=point_new(1,:)-center_ell(1);
point_back(2,:)=point_new(2,:)-center_ell(2);
point_back(3,:)=point_new(3,:)-center_ell(3);

point_back=(evecs_ell')*point_back;

[ ~, alt_back ] = cspice_nearpt( point_back, radii_ell(1), radii_ell(2), radii_ell(3) );

alt_back=reshape(alt_back,s);



