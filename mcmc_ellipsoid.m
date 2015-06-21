ccc

%%
% Next, create a data structure for the observations and control
% variables. Typically one could make a structure |data| that
% contains fields |xdata| and |ydata|.

filename='../CeresShapeModel/SPG/2015_02_26.Ceres-RC2.SPG.2ppd.gaps.xyz';
data=load(filename);

[lon,lat,rad]=cart2sph(data(:,1),data(:,2),data(:,3));

Npts=size(lon,1);

Npts_red=10000;

ind=fix(rand(1,Npts_red)*Npts)+1;

data.lon = lon(ind); 
data.lat = lat(ind); 
data.rad = rad(ind);

a0=500000;
b0=500000;
c0=450000;

param=[a0 b0 c0];

rad=TriEllRadVec(lat,lon,a0,b0,c0,'rad');

cost_fun = @(radii,data) sum((TriEllRadVec(data.lat,data.lon,...
                                    radii(1),radii(2),radii(3),'rad')-data.rad).^2);
                                      
step_param = @(radii,s) [radii(1)+s*(rand-0.5) radii(2)+s*(rand-0.5) radii(3)+s*(rand-0.5)];

N=100000;
s=100;


cost_fun_h = @(param,data) cost_fun(param,data);
step_param_h = @(param,s) step_param(param,s);

mcmc_test(N, s, data, cost_fun_h,step_param_h)

figure; hold on;
hist(param(:,1)/1000,20);

figure; hold on;
hist(param(:,2)/1000,20);

figure; hold on;
hist(param(:,3)/1000,20);