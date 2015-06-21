
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
%  FileName='VestaHASTALAVESTAshape_geoc_elev_3deg_mean_g.grd';

[ri_shape,fii,lambdai]=ReadRawGrid(FileName,{'lon', 'lat', 'z'});

s=size(ri_shape);

[x,y,z]=sph2cart(lambdai,fii,ri_shape*1000);


a=280871;
c=226505; %%  ellipsoid

e=Eccentricity(a,c);
[B,L,H]=XYZ2BLH(x,y,z,a,e);


clear lambdai fii ri_shape

figure('color','k','Position',[1 1 2560 1440]);

Albedo=ones(s(1),s(2));


% Tr1=1;
% Tr2=100;
% 
% 
% h=surf(x(Tr1:Tr2,:),y(Tr1:Tr2,:),z(Tr1:Tr2,:),Albedo(Tr1:Tr2,:));

h=surf(x,y,z,Albedo);
r=sqrt(x.*x+y.*y+z.*z);
% Nstars=1000000;
% 
% [fi_star,lambda_star]=GenerateRandomSphericalCoord(Nstars);
% 
% [x_star,y_star,z_star]=sph2cart(lambda_star,fi_star,1e7);

lighting phong
shading interp
axis equal tight off
box off
hold on;
light_handle=light('Style','infinite','Position',[1 0 0]);

set(h,'BackFaceLighting','lit');
colormap jet;
material([0 0.9 0]);
set(gca, 'Position', [0 0 1 1]);



% plot_star=plot3(x_star,y_star,z_star,'.w','markersize',1);
% 
% r_star=[x_star,y_star,z_star];


% r=1500000;
% axis([-r r -r r -r r]);


set(gca, 'Position', [0 0 1 1]);
% T = get(gca,'tightinset');
% set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

view(0,-75)

%% Reflectance
% ReflFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/DLR/120427_DLR-HAMO-mosp-equi180.png';

% ReflFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/DLR/global_color.png';
% [A,map,alpha]=imread(ReflFileName,'png');
% refl=double(A);
%  
% sr=size(refl);
% 
% lambda=linspace(-180,180,sr(2));
% fi=linspace(-90,90,sr(1));
% 
% [lambdai,fii]=meshgrid(lambda,fi);
%  
% WriteXYZ(lambdai,fii,refl(:,:,1),'global_color_1.xyz');
% WriteXYZ(lambdai,fii,refl(:,:,2),'global_color_2.xyz');
% WriteXYZ(lambdai,fii,refl(:,:,3),'global_color_3.xyz');



% [lambdai,fii,refl]=ReadGRD('/Users/antonermakov/Dawn/Vesta_Shape_Models/DLR/VestaRefl_3m.grd');

[~,~,Rm]=ReadGRD('global_color_1_3m.grd');
[~,~,Gm]=ReadGRD('global_color_2_3m.grd');
[~,~,Bm]=ReadGRD('global_color_3_3m.grd');

Rm=flipud(Rm)';
Gm=flipud(Gm)';
Bm=flipud(Bm)';

Rm=circshift(Rm,[0 3600]);
Gm=circshift(Gm,[0 3600]);
Bm=circshift(Bm,[0 3600]);

RGB=Rm;
RGB(:,:,2)=Gm;
RGB(:,:,3)=Bm;



% refl=refl/255;


% figure;
% pcolor(r(1:10:end,1:10:end)); shading interp;
% refl2=circshift(refl,[0 3600]);
% figure;
% pcolor(flipud(refl2(1:10:end,1:10:end))'); shading interp;
% refl3=flipud(refl2)';
% refl3(refl3==0)=0.52;
% set(h,'CData',refl3);

% refl2=circshift(refl,[0 3600]);

figure;
pcolor((RGB(1:10:end,1:10:end,1))'); shading interp;

set(h,'CData',RGB);


% Reflectance=imread(ReflFileName);
% 
% Reflectance=((double(Reflectance)-1)*0.6/254)';
% 
% s=size(Reflectance);
% 
% fi=linspace(-90,90,s(1));
% lambda=linspace(0,360,s(2));
% [lambdai,fii]=meshgrid(lambda,fi);
% 
% WriteXYZ(lambdai,fii,Reflectance,'Reflectance.xyz');



% campos([-1.366e5, -2.103e5, -1.39e5]*1.1)
% camtarget([0,0,-200000])
% camproj('perspective')
% 
% camup([0 0 -1]);


% J=histeq(refl3);
% 
% set(h,'CData',J);


