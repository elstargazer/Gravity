% FileName='~/Dawn/Balmino2/VestaTest/SH_VestaHASTALAVESTAshape_6min';
FileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POST_VESTA/SHAPE.TXT';

data=load(FileName);

x=data(:,1)*1000;
y=data(:,2)*1000;
z=data(:,3)*1000;

[lambda,fi,r]=cart2sph(x,y,z);

a_approx=280000;
b_approx=280000;
c_approx=226000;

ell2_approx=[a_approx c_approx];
ell3_approx=[a_approx b_approx c_approx];

[ell2_estim,rms2]=FindRotationalEllipsoidRaw(r,fi,lambda,ell2_approx);


% [ell2_estim_spice,rms2_spice]=FindRotationalEllipsoidSpiceRaw(r,fi,lambda,ell2_approx);

% [ell3_estim_spice,rms3_spice]=FindTriaxialEllipsoidSpiceRaw(r,fi,lambda,ell3_approx);


[ell3_estim_spice,rms3_spice]=FindTriaxialEllipsoidSpice(r,fi,lambda,ell3_approx);
[ell2_estim_spice,rms2_spice]=FindRotationalEllipsoidSpice(r,fi,lambda,ell2_approx);


% lmcosi_shape=ReadBalminoSH2(FileName);
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
% 
% 
% a_approx=280000;
% b_approx=280000;
% c_approx=226000;
% 
% ell=[a_approx b_approx c_approx];
% 
% [ell_estim,rms]=FindTriaxialEllipsoidSpice(ri_shape,fii,lambdai,ell);