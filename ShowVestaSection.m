
%% Load spherical harmonics shape model
% load VestaHASTALAVESTAshape_sh720.mat
% 
% lmcosi_shape=ReadBalminoSH2('~/Dawn/Balmino2/VestaTest/SH_VestaHASTALAVESTAshape_6min');
% 
% Resolution=1;
% MaxDegreeTopo=100;
% 
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1); 
%  
% [ri_shape,~]=plm2xyz(lmcosi_shape,Resolution);
% 
% s=size(ri_shape);
% 
% fi=linspace(90,-90,s(1))/180*pi;
% lambda=linspace(0,360,s(2))/180*pi;
% 
% [lambdai,fii]=meshgrid(lambda,fi);
% 
% [x,y,z]=sph2cart(lambdai,fii,ri_shape);

%% Load grid shape model

FileName='VestaHASTALAVESTAshape_geoc_elev_3min_mean_g.grd';
% FileName='VestaHASTALAVESTAshape_geoc_elev_3deg_mean_g.grd';

[ri_shape,fii,lambdai]=ReadRawGrid(FileName);


s=size(ri_shape);

[x,y,z]=sph2cart(lambdai,fii,ri_shape*1000);



% Tr1=1;
% Tr2=Tr1+s(2)-1;


Tr1=1;
Tr2=fix(0.7438*s(1));

% Tr1=Tr2;
% Tr2=s(1);

lambda1=lambdai(Tr1,1)*180/pi;
lambda2=lambdai(Tr2,1)*180/pi;

clear lambdai fii ri_shape

figure('color','k','Position',[1 1 2560 1440]);

Albedo=ones(s(1),s(2));


h=surf(x(Tr1:Tr2,:),y(Tr1:Tr2,:),z(Tr1:Tr2,:),Albedo(Tr1:Tr2,:));

% h=surf(x,y,z,Albedo);

lighting phong
shading interp
axis equal tight off
box off

hold on;
light_handle=light('Style','infinite');

set(h,'BackFaceLighting','lit');

colormap gray;

material([0 0.9 0]);


set(gca, 'Position', [0 0 1 1]);
% T = get(gca,'tightinset');
% set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

% view(240,10)
view(210,-15)


fi_f1=(-90:1:90)';
fi_f2=(90:-1:-90)';
lambda_f1=ones(numel(fi_f1),1).*lambda1;
lambda_f2=ones(numel(fi_f2),1).*lambda2;

fi_f=[fi_f1];
lambda_f=[lambda_f1];


r_shape_f=plm2xyz(lmcosi_shape,fi_f(:),lambda_f(:));
r_mantle_f=plm2xyz(lmcosi_mantle_shape_new,fi_f(:),lambda_f(:));

[x_shape_f,y_shape_f,z_shape_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,r_shape_f);

[x_mantle_f,y_mantle_f,z_mantle_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,r_mantle_f);

lambda_f_perp=(lambda1)-90;
[x_perp,y_perp,z_perp]=sph2cart(lambda_f_perp/180*pi,0,1000);

patch_crust=patch(x_shape_f,y_shape_f,z_shape_f,[0.9 0.9 0.9],'AmbientStrength',0.4,'SpecularStrength',1,'DiffuseStrength',1);
patch_mantle=patch(x_mantle_f+x_perp,y_mantle_f+y_perp,z_mantle_f+z_perp,[0.5 0.7 0.5],'AmbientStrength',0.4,'SpecularStrength',0,'DiffuseStrength',1)%,'EdgeColor','none');

fi_light=0;
lambda_light=-40;

[x_light,y_light,z_light]=sph2cart(lambda_light/180*pi,fi_light/180*pi,1);


set(light_handle,'Position',[x_light y_light z_light]);

[x_core_f,y_core_f,z_core_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,1);

x_core_f=x_core_f.*a_core;
y_core_f=y_core_f.*b_core;
z_core_f=z_core_f.*c_core;

patch_core=patch(x_core_f+2*x_perp,y_core_f+2*y_perp,z_core_f+2*z_perp,[.95 .95 .95],'AmbientStrength',0.4,'SpecularStrength',0.5,'DiffuseStrength',0.1)%,'EdgeColor','none');

fi_f=[fi_f2];
lambda_f=[lambda_f2];


r_shape_f=plm2xyz(lmcosi_shape,fi_f(:),lambda_f(:));
r_mantle_f=plm2xyz(lmcosi_mantle_shape_new,fi_f(:),lambda_f(:));

[x_shape_f,y_shape_f,z_shape_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,r_shape_f);

[x_mantle_f,y_mantle_f,z_mantle_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,r_mantle_f);

lambda_f_perp=(lambda2)-270;
[x_perp,y_perp,z_perp]=sph2cart(lambda_f_perp/180*pi,0,1000);

patch_crust2=patch(x_shape_f,y_shape_f,z_shape_f,[0.9 0.9 0.9],'AmbientStrength',0.4,'SpecularStrength',1,'DiffuseStrength',1);
patch_mantle2=patch(x_mantle_f+x_perp,y_mantle_f+y_perp,z_mantle_f+z_perp,[0.5 0.7 0.5],'AmbientStrength',0.4,'SpecularStrength',0,'DiffuseStrength',1)%,'EdgeColor','none');


[x_core_f,y_core_f,z_core_f]=sph2cart(lambda_f/180*pi,fi_f/180*pi,1);

x_core_f=x_core_f.*a_core;
y_core_f=y_core_f.*b_core;
z_core_f=z_core_f.*c_core;

patch_core2=patch(x_core_f+2*x_perp,y_core_f+2*y_perp,z_core_f+2*z_perp,[.95 .95 .95],'AmbientStrength',0.4,'SpecularStrength',0.5,'DiffuseStrength',0.1)%,'EdgeColor','none');


