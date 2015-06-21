
%% Load shaoe model in SH
% MaxDegreeTopo=720;
% 
% load VestaHASTALAVESTAshape_sh1500.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
% % Create a mesh
% Step=.1;
% [lambdai,fii]=meshgrid(0:Step:360,90:-Step:-90);
% lambdai=lambdai/180*pi;
% fii=fii/180*pi;
% ri=plm2xyz(lmcosi_shape,Step);
% 
% 
% [xi,yi,zi]=sph2cart(lambdai,fii,ri);

%% Or load grided shape model

ShapeModelFileName1='~/Dawn/Balmino/VestaTest/VestaPOSTVESTAshape_geoc_elev12m.grd';
ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR12m.grd';

[lambdai,fii,ri]=ReadGRD(ShapeModelFileName2);

% ri=ri*1000;
lambdai=lambdai/180*pi;
fii=fii/180*pi;

[xi,yi,zi]=sph2cart(lambdai,fii,ri);



xi=xi/1000;
yi=yi/1000;
zi=zi/1000;

%% Computing Curvature

[gm, samc] = mcurvature_vec(xi,yi,zi);

% shift=1801;
% 
% xi2=circshift(xi,[0 shift]);
% yi2=circshift(yi,[0 shift]);
% zi2=circshift(zi,[0 shift]);

% [gm, samc] = mcurvature_vec(xi,yi2,zi2);
% 
% gm3=circshift(gm2,[0 -shift]);
figure; hold on;
R=1./gm;
log10absR=log10(abs(R));
hist(log10absR(:),100);


% MapRadialGrid(R,-500,500);

Minlog10absR=2.5;
Maxlog10absR=7;

% MapRadialGrid(log10absR,Minlog10absR,Maxlog10absR);

log10absR(log10absR<Minlog10absR)=2.0;

% CurvGridFileName='~/Dawn/Balmino/VestaTest/VestaCurvGrid.xyz';
% WriteXYZ(lambdai*180/pi,fii*180/pi,log10absR,CurvGridFileName);


%% Curv in SH

[lmcosi_R]=xyz2plm(log10absR,20,'im',[],[],[]);

[log10absR_sh,lon,lat]=plm2xyz(lmcosi_R,12/60);

[lon,lat]=meshgrid(lon,lat);

MapRadialGrid(log10absR_sh,Minlog10absR,Maxlog10absR);
caxis([2 2.1]);


figure; hold on;
hist(log10absR_sh(:),50);
 
% CurvGridFileName_sh='~/Dawn/Balmino/VestaTest/VestaCurvGrid_sh.xyz';
% WriteXYZ(lon,lat,log10absR_sh,CurvGridFileName_sh);



