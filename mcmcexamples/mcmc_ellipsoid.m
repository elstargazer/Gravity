ccc

%%
% Next, create a data structure for the observations and control
% variables. Typically one could make a structure |data| that
% contains fields |xdata| and |ydata|.

filename='~/Dawn/CeresShapeModel/SPG/2015_02_26.Ceres-RC2.SPG.2ppd.gaps.xyz';
data=load(filename);

[lon,lat,rad]=cart2sph(data(:,1),data(:,2),data(:,3));

[x,y,z]=sph2cart(lon,lat,rad);

[ center, radii, evecs, v ] = ellipsoid_fit( [x y z], 0 );

x=x-center(1);
y=y-center(2);
z=z-center(3);

[ center, radii, evecs, v ] = ellipsoid_fit( [x y z], 1 );

Npts=size(lon,1);
Npts_red=1000;

ind=fix(rand(1,Npts_red)*Npts)+1;

data.lon = lon(:); 
data.lat = lat(:); 
data.rad = rad(:);

a0=radii(1)+randn*1000;
b0=radii(2)+randn*1000;
c0=radii(3)+randn*1000;

rad=TriEllRadVec(lat,lon,a0,b0,c0,'rad');

cost_fun = @(radii,data) sum((TriEllRadVec(data.lat,data.lon,...
                              radii(1),radii(2),radii(3),'rad')-data.rad).^2);
                               
step_param = @(radii,s) [radii(1)+s*(rand-0.5) radii(2)+s*(rand-0.5) radii(3)+s*(rand-0.5)];

N=100000;
s=100;

cost_fun_h = @(param,data) cost_fun(param,data);
step_param_h = @(param,s) step_param(param,s);

param0=[a0 b0 c0];

Nc=8;

param=cell(1,Nc);

tic
parfor j=1:Nc
    param{j}=mcmc_test(N, s, data, param0, cost_fun_h,step_param_h);    
end
toc

param_all=[];

for i=1:Nc    
    param_all=[param_all; param{i}];    
end
    
figure;

subplot(3,1,1); hold on;
hist(param_all(:,1)/1000,20);
title('a','FontSize',12);
am=mean(param_all(:,1))/1000;

subplot(3,1,2); hold on;
hist(param_all(:,2)/1000,20);
title('b','FontSize',12);
bm=mean(param_all(:,2))/1000;

subplot(3,1,3); hold on;
hist(param_all(:,3)/1000,20);
title('c','FontSize',12);
cm=mean(param_all(:,3))/1000;

[am bm cm]

radii'/1000
