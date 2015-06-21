cc

a=278000;
b=284000;
c=226000;

radii_ell_ini=[a b c];

[fi,lambda]=meshgrid(-pi/2:pi/100:pi/2,0:pi/100:2*pi);


x_mesh=a*cos(fi).*cos(lambda);
y_mesh=b*cos(fi).*sin(lambda);
z_mesh=c*sin(fi);

% surf(x,y,z);
% axis equal;
% hold on;


s=size(x_mesh);

sigma=10000;



x=x_mesh(:)+sigma*randn(size(x_mesh(:)));
y=y_mesh(:)+sigma*randn(size(y_mesh(:)));
z=z_mesh(:)+sigma*randn(size(z_mesh(:)));



% x=x_mesh(:)+sigma*sin(3*fi(:));
% y=y_mesh(:)+sigma*sin(3*fi(:));
% z=z_mesh(:)+sigma*sin(3*fi(:));


% x_mesh=reshape(x,s);
% y_mesh=reshape(y,s);
% z_mesh=reshape(z,s);



% surf(x_mesh,y_mesh,z_mesh);
% axis equal

point_ini=[x, y, z]';



% surf(x,y,z);
% axis equal


[ center_ell, radii_ell, evecs_ell, v_ell ] = ellipsoid_fit( point_ini', 0 );


point=(evecs_ell')*point_ini;

[ npoint, alt ] = cspice_nearpt( point, radii_ell(1), radii_ell(2), radii_ell(3) );


%% Transformation

R=rot(pi/4,1)*rot(pi/6,2)*rot(pi/3,3);

point_new=R*point;

x_new=point_new(1,:);
y_new=point_new(2,:);
z_new=point_new(3,:);

[ center_ell_new, radii_ell_new, evecs_ell_new, v_ell_new ] = ellipsoid_fit( point_new', 0 );


%% Transformation back

% evecs_ell_new2=[-evecs_ell_new(:,2) evecs_ell_new(:,1) evecs_ell_new(:,3)];
% evecs_ell_new=evecs_ell_new2;


%  point_back=(R')*point_new;
 point_back=(evecs_ell_new')*point_new;

[ center_ell_back, radii_ell_back, evecs_ell_back, v_ell_back ] = ellipsoid_fit( point_back', 0 );
[ npoint_back, alt_back ] = cspice_nearpt( point_back, radii_ell_back(1), radii_ell_back(2), radii_ell_back(3) );
% [ npoint_back, alt_back ] = cspice_nearpt( point_back, a, b, c );
% [ npoint_back, alt_back ] = cspice_nearpt( point_back, a, b, c );

x_back=point_back(1,:);
y_back=point_back(2,:);
z_back=point_back(3,:);


x_back=reshape(x_back,s);
y_back=reshape(y_back,s);
z_back=reshape(z_back,s);

% surf(x_back,y_back,z_back);

figure
plot(abs(alt-alt_back),'.');
 
logdiff=log10(abs((alt-alt_back)./alt));
 
% hist(logdiff,40);
% plot(logdiff)

%  plot(alt,'.')






