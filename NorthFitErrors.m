ccc;

%% Load shape model
ShapeModelFile='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POST_VESTA/SHAPE.TXT';

data=load(ShapeModelFile);

x_raw=data(:,1);
y_raw=data(:,2);
z_raw=data(:,3);

N_raw=numel(x_raw);

[lambda_raw,fi_raw,r_raw]=cart2sph(x_raw,y_raw,z_raw);


MaxDegreeTopo=100;
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

[ri_shape,lambda,fi,~]=plm2xyz(lmcosi_shape,1);

[lambdai,fii]=meshgrid(lambda/180*pi,fi/180*pi);


[x,y,z]=sph2cart(lambdai,fii,ri_shape);

[~,~,H]=XYZ2BLH(x,y,z,280917.27,Eccentricity(280917.27,226248.16));


% vesta_fig=figure;
% surf(x/1000,y/1000,z/1000,H);
% shading interp
% axis equal
% light
% lighting phong
% axis tight equal off

% 

% MapRadialGrid(flipud(H)/1000);
% hold on
open CrustalThicknessFigure.fig;

%% Draw boundary

NA=10

for Try=1:NA

    Try

Boundary=inputm();

% load BoundaryNorthFit.mat

Boundary(end+1,:)=Boundary(1,:);
% [lattrk,lontrk] = track(fliplr(Boundary)*180/pi);

[lattrk,lontrk] = track((Boundary)*180/pi);

WriteXYZ(lontrk,lattrk,[lontrk*0+1],'North_Boundary.xyz');

lattrk(isnan(lattrk))=[];
lontrk(isnan(lontrk))=[];

%  plotm(lattrk/180*pi,lontrk/180*pi,'k-','LineWidth',5);


fi_raw_bound=BoundaryLatitude(lambda_raw*180/pi,lattrk,lontrk);

plotm(lattrk/180*pi,lontrk/180*pi,'-k','LineWidth',1);


Condition=(fi_raw>fi_raw_bound/180*pi);
x_raw_cut=x_raw(Condition);
y_raw_cut=y_raw(Condition);
z_raw_cut=z_raw(Condition);

N_raw_cut=numel(x_raw_cut);


[lambda_raw_cut,fi_raw_cut,r_raw_cut]=cart2sph(x_raw_cut,y_raw_cut,z_raw_cut);

% plotm(fi_raw_cut,lambda_raw_cut,'.k','markersize',1);


N_rand_pts=100000;

index_rand=fix(rand(1,N_rand_pts)*N_raw);

index_rand=index_rand(index_rand>0);

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

center_cuti(:,Try)=center_cut;
radii_cuti(:,Try)=radii_cut;
evecs_cuti(:,:,Try)=evecs_cut;

a_vec_cuti(:,Try)=a_vec_cut;
b_vec_cuti(:,Try)=b_vec_cut;
c_vec_cuti(:,Try)=c_vec_cut;


end

dx_cuti=center_cuti(1,:);
dy_cuti=center_cuti(2,:);
dz_cuti=center_cuti(3,:);


a_cuti=radii_cuti(1,:);
b_cuti=radii_cuti(2,:);
c_cuti=radii_cuti(3,:);


lambda_a_cuti=atan2(a_vec_cuti(2,:),a_vec_cuti(1,:))*180/pi;
lambda_b_cuti=atan2(b_vec_cuti(2,:),b_vec_cuti(1,:))*180/pi;
lambda_c_cuti=atan2(c_vec_cuti(2,:),c_vec_cuti(1,:))*180/pi;


fi_a_cuti=atan2(a_vec_cuti(3,:),sqrt((a_vec_cuti(1,:)).^2+(a_vec_cuti(2,:)).^2))*180/pi;
fi_b_cuti=atan2(b_vec_cuti(3,:),sqrt((b_vec_cuti(1,:)).^2+(b_vec_cuti(2,:)).^2))*180/pi;
fi_c_cuti=atan2(c_vec_cuti(3,:),sqrt((c_vec_cuti(1,:)).^2+(c_vec_cuti(2,:)).^2))*180/pi;


%% Mean values
% 
% lambda_a_cuti=To0360Range(lambda_a_cuti);
% lambda_b_cuti=To0360Range(lambda_b_cuti);
% lambda_c_cuti=To0360Range(lambda_c_cuti);


% fi_a_cuti=To0360Range(fi_a_cuti);
% fi_b_cuti=To0360Range(fi_b_cuti);
% fi_c_cuti=To0360Range(fi_c_cuti);



a_cut_mean=mean(a_cuti);
b_cut_mean=mean(b_cuti);
c_cut_mean=mean(c_cuti);

a_eq_cut=sqrt(a_cut_mean*b_cut_mean);
f_cut=(a_eq_cut-c_cut_mean)/a_eq_cut

figure;hold on; plot(a_cuti,'-r'); plot(b_cuti,'-g'); plot(c_cuti,'-b');

axes_mean=[a_cut_mean b_cut_mean c_cut_mean]

dx_cut_mean=mean(dx_cuti);
dy_cut_mean=mean(dy_cuti);
dz_cut_mean=mean(dz_cuti);


figure;hold on; plot(dx_cuti,'-r'); plot(dy_cuti,'-g'); plot(dz_cuti,'-b');

dr_mean=[dx_cut_mean dy_cut_mean dz_cut_mean]

lambda_a_cut_mean=mean(lambda_a_cuti);
lambda_b_cut_mean=mean(lambda_b_cuti);
lambda_c_cut_mean=mean(lambda_c_cuti);

figure;hold on; plot(lambda_a_cuti,'-r'); plot(lambda_b_cuti,'-g'); plot(lambda_c_cuti,'-b');

lambda_abc_mean=[lambda_a_cut_mean lambda_b_cut_mean lambda_c_cut_mean]

fi_a_cut_mean=mean(fi_a_cuti);
fi_b_cut_mean=mean(fi_b_cuti);
fi_c_cut_mean=mean(fi_c_cuti);

figure;hold on; plot(fi_a_cuti,'-r'); plot(fi_b_cuti,'-g'); plot(fi_c_cuti,'-b');


fi_abc_mean=[fi_a_cut_mean fi_b_cut_mean fi_c_cut_mean]

%% Uncertanties

a_cut_std=std(a_cuti);
b_cut_std=std(b_cuti);
c_cut_std=std(c_cuti);

abc_std=[a_cut_std b_cut_std c_cut_std]

dx_cut_std=std(dx_cuti);
dy_cut_std=std(dy_cuti);
dz_cut_std=std(dz_cuti);

dr_std=[dx_cut_std dy_cut_std dz_cut_std]

lambda_a_cut_std=std(lambda_a_cuti);
lambda_b_cut_std=std(lambda_b_cuti);
lambda_c_cut_std=std(lambda_c_cuti);

lambda_abc_std=[lambda_a_cut_std lambda_b_cut_std lambda_c_cut_std]


fi_a_cut_std=std(fi_a_cuti);
fi_b_cut_std=std(fi_b_cuti);
fi_c_cut_std=std(fi_c_cuti);

fi_abc_std=[fi_a_cut_std fi_b_cut_std fi_c_cut_std]












































