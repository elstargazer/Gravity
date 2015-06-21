ccc

%% Parameters
MaxDegree=2;

ag=290000;
cg=265000;

Step=1;


%% Load gravity model
FileName='/Users/antonermakov/Dawn/Gravity/Frank20/grvfld.dawn.img.test37';

[lmcosi,mu,Rref]=ReadGEODYNGravity(FileName);
lmcosi=TruncateGravityModel(lmcosi,MaxDegree,1);

%% Choose reference surface (ellipsoid of revolition

lambda=(0:Step:360)/180*pi;
fi=(-90:Step:90)/180*pi;

[fii,lambdai]=meshgrid(fi,lambda);
[xr,yr,zr]=TwoEllRadVec(lambdai,fii,ag,cg);

%% Compute gravity at the reference surface

[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi,xr,yr,zr);
[g_up_obs,g_east_obs,g_north_obs]=GravityComponents(gx_obs,gy_obs,gz_obs,xr,yr,zr,ag,cg);

%% Plotting results

AGUaxes;
pcolorm(fii,lambdai,g_up_obs*1e5); shading interp;
cbar=colorbar('FontSize',20);
ylabel(cbar,'Gravity Acceleration [mGal]','FontSize',20);



