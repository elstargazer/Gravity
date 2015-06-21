cc
%% Load shape model
ShapeModelFile='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POST_VESTA/SHAPE.TXT';

data=load(ShapeModelFile);

x_raw=data(:,1);
y_raw=data(:,2);
z_raw=data(:,3);

N_raw=numel(x_raw);

[lambda_raw,fi_raw,r_raw]=cart2sph(x_raw,y_raw,z_raw);


MaxDegreeTopo=500;
Resolution=0.25
load VestaHASTALAVESTAshape_sh720.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,Resolution);

[ri_shape,lambda,fi,~]=plm2xyz(lmcosi_shape,Resolution);

[lambdai,fii]=meshgrid(lambda/180*pi,fi/180*pi);


[x,y,z]=sph2cart(lambdai,fii,ri_shape);

[~,~,H]=XYZ2BLH(x,y,z,280917.27,Eccentricity(280917.27,226248.16));


vesta_fig=figure;
surf(x/1000,y/1000,z/1000,H);
shading interp
axis equal
light
lighting phong
axis tight equal off

% figure
% 
% MapRadialGrid(flipud(H)/1000);
% hold on

open CrustalThicknessFigure.fig;

%% Draw boundary

Boundary=inputm();

% load BoundaryNorthFit.mat

Boundary(end+1,:)=Boundary(1,:);
% [lattrk,lontrk] = track(fliplr(Boundary)*180/pi);

[lattrk,lontrk] = track((Boundary)*180/pi);

plotm(lattrk/180*pi,lontrk/180*pi,'-k','LineWidth',3);

lattrk(isnan(lattrk))=[];
lontrk(isnan(lontrk))=[];

WriteXYZ(lontrk,lattrk,[lontrk*0+1],'North_Boundary.xyz');


% plotm(lattrk/180*pi,lontrk/180*pi,'k-','LineWidth',5);


fi_raw_bound=BoundaryLatitude(lambda_raw*180/pi,lattrk,lontrk);

Condition=(fi_raw>fi_raw_bound/180*pi);
x_raw_cut=x_raw(Condition);
y_raw_cut=y_raw(Condition);
z_raw_cut=z_raw(Condition);

N_raw_cut=numel(x_raw_cut);


[lambda_raw_cut,fi_raw_cut,r_raw_cut]=cart2sph(x_raw_cut,y_raw_cut,z_raw_cut);

% plotm(fi_raw_cut,lambda_raw_cut,'.k','markersize',1);


N_rand_pts=10000;

index_rand=fix(rand(1,N_rand_pts)*N_raw);

x_raw_rand_sel=x_raw(index_rand);
y_raw_rand_sel=y_raw(index_rand);
z_raw_rand_sel=z_raw(index_rand);


%  [ center, radii, evecs, v ] = ellipsoid_fit( [x_raw_rand_sel y_raw_rand_sel z_raw_rand_sel], 0 );
 
[ center, radii, evecs, v ] = ellipsoid_fit( [x_raw y_raw z_raw], 0 );
 
radii_rot=[sqrt(radii(1)*radii(2)) sqrt(radii(1)*radii(2)) radii(3)];

f=(radii_rot(1)-radii_rot(3))/radii_rot(1);
 
a_vec=evecs(:,1);
b_vec=evecs(:,2);
c_vec=evecs(:,3);
 
[ center_cut, radii_cut, evecs_cut, v_cut ] = ellipsoid_fit( [x_raw_cut y_raw_cut z_raw_cut], 0 );
 
% save('NorthEllipsoidOrientation.mat','center_cut', 'radii_cut', 'evecs_cut');


radii_rot_cut=[sqrt(radii_cut(1)*radii_cut(2)) sqrt(radii_cut(1)*radii_cut(2)) radii_cut(3)];

f_cut=(radii_rot_cut(1)-radii_rot_cut(3))/radii_rot_cut(1);
 
a_vec_cut=evecs_cut(:,1);
b_vec_cut=evecs_cut(:,2);
c_vec_cut=evecs_cut(:,3);
  
 
figure(vesta_fig)
hold on

%% Draw fit
maxd = max( radii );
step = maxd / 50;


[ x_fit, y_fit, z_fit ] = meshgrid( -maxd:step:maxd + center_cut(1), -maxd:step:maxd + center_cut(2), -maxd:step:maxd + center_cut(3) );

Ellipsoid = v_cut(1) *x_fit.*x_fit +   v_cut(2) * y_fit.*y_fit + v_cut(3) * z_fit.*z_fit + ...
          2*v_cut(4) *x_fit.*y_fit + 2*v_cut(5)*x_fit.*z_fit + 2*v_cut(6) * y_fit.*z_fit + ...
          2*v_cut(7) *x_fit    + 2*v_cut(8)*y_fit    + 2*v_cut(9) * z_fit;
      
p = patch( isosurface( x_fit, y_fit, z_fit, Ellipsoid, 1 ) );
set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
view( -70, 40 );
axis vis3d;


view(0,0)

zlim([-300.000 300.000]);
xlim([-300.000 300.000]);
ylim([-300 300]);
 

% plotplm(lmcosi_shape,[],[],2,1,[],[],[]);    

hold on

vec_lenght=500000;
 
%% Ellipsoid axis for the North shape
line([center_cut(1) center_cut(1)+vec_lenght*a_vec_cut(1)],...
    [center_cut(2) center_cut(2)+vec_lenght*a_vec_cut(2)],...
    [center_cut(3) center_cut(3)+vec_lenght*a_vec_cut(3)],'Color','r','LineWidth',3,'LineStyle','-');

% line([center_cut(1) center_cut(1)+vec_lenght*b_vec_cut(1)],...
%     [center_cut(2) center_cut(2)+vec_lenght*b_vec_cut(2)],...
%     [center_cut(3) center_cut(3)+vec_lenght*b_vec_cut(3)],'Color','g','LineWidth',3,'LineStyle','-');
% 
% line([center_cut(1) center_cut(1)+vec_lenght*c_vec_cut(1)],...
%     [center_cut(2) center_cut(2)+vec_lenght*c_vec_cut(2)],...
%     [center_cut(3) center_cut(3)+vec_lenght*c_vec_cut(3)],'Color','b','LineWidth',3,'LineStyle','-');

%% Ellipsoid axis for the Total shape

line([center(1) center(1)+vec_lenght*a_vec(1)],...
    [center(2) center(2)+vec_lenght*a_vec(2)],...
    [center(3) center(3)+vec_lenght*a_vec(3)],'Color','r','LineWidth',3,'LineStyle','--');

% line([center(1) center(1)+vec_lenght*b_vec(1)],...
%     [center(2) center(2)+vec_lenght*b_vec(2)],...
%     [center(3) center(3)+vec_lenght*b_vec(3)],'Color','g','LineWidth',3,'LineStyle','--');
% 
% line([center(1) center(1)+vec_lenght*c_vec(1)],...
%     [center(2) center(2)+vec_lenght*c_vec(2)],...
%     [center(3) center(3)+vec_lenght*c_vec(3)],'Color','b','LineWidth',3,'LineStyle','--');

%% Rotation axis

line([0 0],...
    [0 0],...
    [-vec_lenght vec_lenght],'Color','k','LineWidth',4,'LineStyle','-');

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

alpha(0.7)

f_cut

lambda_a=atan2(a_vec(2),a_vec(1))*180/pi
lambda_b=atan2(b_vec(2),b_vec(1))*180/pi
lambda_c=atan2(c_vec(2),c_vec(1))*180/pi

lambda_a_cut=atan2(a_vec_cut(2),a_vec_cut(1))*180/pi
lambda_b_cut=atan2(b_vec_cut(2),b_vec_cut(1))*180/pi
lambda_c_cut=atan2(c_vec_cut(2),c_vec_cut(1))*180/pi

fi_a=atan2(a_vec(3),norm([a_vec(1) a_vec(2)]))*180/pi
fi_b=atan2(b_vec(3),norm([b_vec(1) b_vec(2)]))*180/pi
fi_c=atan2(c_vec(3),norm([c_vec(1) c_vec(2)]))*180/pi



fi_a_cut=atan2(a_vec_cut(3),norm([a_vec_cut(1) a_vec_cut(2)]))*180/pi
fi_b_cut=atan2(b_vec_cut(3),norm([b_vec_cut(1) b_vec_cut(2)]))*180/pi
fi_c_cut=atan2(c_vec_cut(3),norm([c_vec_cut(1) c_vec_cut(2)]))*180/pi
