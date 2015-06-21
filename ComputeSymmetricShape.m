
Resolution=0.1;
MaxDegreeTopo=500;
 
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

[ri_shape,~]=plm2xyz(lmcosi_shape,Resolution);

s=size(ri_shape);

fi=linspace(90,-90,s(1))/180*pi;
lambda=linspace(0,360,s(2))/180*pi;


[lambdai,fii]=meshgrid(lambda,fi);

ell=[285000 22600];

% [ell,ell_rms]=FindRotationalEllipsoid(ri_shape,fii,lambdai,ell);

% [x,y,z]=sph2cart(lambdai,fii,ri_shape);

% figure
% surf(x,y,z,ri_shape)
% axis equal
% f=(ell(1)-ell(2))/ell(1);




half=(s(1)-1)/2;


%% reflect North to South

ri_shape(half+2:end,:)=flipud(ri_shape(1:half,:));

[x,y,z]=sph2cart(lambdai,fii,ri_shape);

figure('Color',[0 0 0]);
surf(x,y,z,ri_shape)
StandardLight

% [ell_refl,ell_rms_refl]=FindRotationalEllipsoid(ri_shape,fii,lambdai,ell);

% f_refl=(ell_refl(1)-ell_refl(2))/ell_refl(1);