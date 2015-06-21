ccc
ag=2.941232378739690e+05;
bg=2.649866853974999e+05;
Rref=265000;

mu=1.728824342710000e+10;



ReadBalminoSHGravity

FiStep=1;
LambdaStep=1;
[fii,lambdai]=meshgrid( -pi/2:FiStep/180*pi:pi/2 , 0:LambdaStep/180*pi:2*pi);



x=ag*cos(lambdai).*cos(fii);
y=ag*sin(lambdai).*cos(fii);
z=bg*sin(fii);

g=RadialGForceFromGravityModel(x,y,z,C_gt,S_gt,mu,Rref)*1e5;

% AGUaxes


% g_step=200; %% in mGals
% 
% levels=fix(min(min(g))/g_step)*g_step:g_step:fix(max(max(g))/g_step)*g_step;
% [C,h]=contourm(fii,lambdai,g,levels,'linewidth',2);
% ht = clabelm(C,h,levels);
% set(ht,'Color','w','BackgroundColor','none','FontWeight','bold','fontsize',16);
% 

GridFileName='gxyz_g_topo_sh_10.txt';
%MapFileName='gmap.ps';


WriteXYZ(lambdai*180/pi,fii*180/pi,g,GridFileName);

