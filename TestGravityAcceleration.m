ccc

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20F/JGV20F02.SHA';

[C_obs,S_obs,C_obs_std,S_obs_std,mu_obs,Rref]=LoadGravityModel(GravityFileName);
[lmcosi_obs,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);




ag=2.918299428885863e+05;
bg=2.650067859489697e+05;


FiStep=1;
LambdaStep=1;
[fii,lambdai]=meshgrid( -pi/2:FiStep/180*pi:pi/2 , 0:LambdaStep/180*pi:2*pi);
x=ag*cos(lambdai).*cos(fii);
y=ag*sin(lambdai).*cos(fii);
z=bg*sin(fii);

r=sqrt(x.*x+y.*y+z.*z);

tic
g_r=1e5*RadialGravityAcceleration(x,y,z,C_obs,S_obs,mu_obs,Rref);
toc
tic
[gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_obs,x,y,z);

gx=gx*1e5;
gy=gy*1e5;
gz=gz*1e5;

toc

g=sqrt(gx.*gx+gy.*gy+gz.*gz);

%% Radial gravity
[fig,ax]=AGUaxes;
MinValue=19000;
MaxValue=25000;
StepLevels=250;
StepLabels=1000;
Levels=MinValue:StepLevels:MaxValue;
[C,h]=contourfm(fii,lambdai,g_r,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')
Labels=MinValue:StepLabels:MaxValue;
ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Radial gravity component','fontsize',25);

%% Total gravity
[fig2,ax2]=AGUaxes;
% MinValue=min(min(g));
% MaxValue=max(max(g));
StepLevels=250;
StepLabels=1000;
Levels=MinValue:StepLevels:MaxValue;
[C,h]=contourfm(fii,lambdai,g,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')
Labels=MinValue:StepLabels:MaxValue;
ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Gravity vector magnitude','fontsize',25);

%% Difference
dg=g-g_r;

[fig3,ax3]=AGUaxes;
StepLevels=25;
StepLabels=50;
MinValue=-mod(min(min(dg)),StepLevels)+min(min(dg));
MaxValue=mod(-max(max(dg)),StepLevels)+max(max(dg));
Levels=MinValue:StepLevels:MaxValue;
Labels=MinValue:StepLabels:MaxValue;
[C,h]=contourm(fii,lambdai,dg,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')

ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Gravity vector magnitude','fontsize',25);


%% Check
g_xu=gx./g;
g_yu=gy./g;
g_zu=gz./g;


xu=x./r;
yu=y./r;
zu=z./r;

l=g_xu.*xu+g_yu.*yu+g_zu.*zu;

g_r2=l.*g;

d_gr=g_r-g_r2;

[fig4,ax4]=AGUaxes;
StepLevels=250;
StepLabels=1000;
MinValue=19000;%-mod(min(min(g_r2)),StepLevels)+min(min(g_r2));
MaxValue=25000;%mod(-max(max(g_r2)),StepLevels)+max(max(g_r2));
Levels=MinValue:StepLevels:MaxValue;
Labels=MinValue:StepLabels:MaxValue;
[C,h]=contourfm(fii,lambdai,g_r2,Levels,'linewidth',5);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')

ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Gravity vector magnitude','fontsize',25);


%% Check 2
dg2=g-g_r2;

[fig5,ax5]=AGUaxes;
StepLevels=25;
StepLabels=50;
MinValue=-mod(min(min(dg2)),StepLevels)+min(min(dg2));
MaxValue=mod(-max(max(dg2)),StepLevels)+max(max(dg2));
Levels=MinValue:StepLevels:MaxValue;
Labels=MinValue:StepLabels:MaxValue;
[C,h]=contourm(fii,lambdai,dg2,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')

ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Gravity vector magnitude','fontsize',25);



%% Ellipsoidal normals

[Nx,Ny,Nz] = surfnorm(x,y,z);


N=sqrt(Nx.*Nx+Ny.*Ny+Nz.*Nz);

Nxu=-Nx./N;
Nyu=-Ny./N;
Nzu=-Nz./N;

l_ell=g_xu.*Nxu+g_yu.*Nyu+g_zu.*Nzu;

g_ell=l_ell.*g;

[fig6,ax6]=AGUaxes;
StepLevels=250;
StepLabels=1000;
MinValue=19000;
MaxValue=25000;
Levels=MinValue:StepLevels:MaxValue;
Labels=MinValue:StepLabels:MaxValue;
[C,h]=contourfm(fii,lambdai,g_ell,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')

ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Gravity vector norma to ellipsoid','fontsize',25);

%% Ellipsoidal - g

dg_ell=g_ell-g;

[fig7,ax7]=AGUaxes;
StepLevels=2.5;
StepLabels=10;
MinValue=-mod(min(min(dg_ell)),StepLevels)+min(min(dg_ell));
MaxValue=mod(-max(max(dg_ell)),StepLevels)+max(max(dg_ell));
Levels=MinValue:StepLevels:MaxValue;
Labels=MinValue:StepLabels:MaxValue;
[C,h]=contourm(fii,lambdai,dg_ell,Levels,'linewidth',2);
% colormap jet;
% axis([axis MinValue MaxValue]);
% colorbar('peer',ax,'fontsize',25,'location','northoutside')

ht = clabelm(C,h,Labels);
set(ht,'Color','k','BackgroundColor','none','FontWeight','bold','fontsize',16);
title('Ell - g ','fontsize',25);

















