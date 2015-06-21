ccc

%% Loading shape models
%% Load Shape model
filename1='/Users/ermakov/Dawn/ShapeModel/Gaskell/20140513/Vesta20140513shape_geoc_elev.grd';
filename2='/Users/ermakov/Dawn/ShapeModel/DLR/DLR_SHAPE_2014_06_19/VestaShapeDLR64.grd';

[lambdai1,fii1,ri1]=ReadGRD(filename1);
[~       ,~   ,ri2]=ReadGRD(filename2); ri2=flipud(ri2/1000);

[xi1,yi1,zi1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);
[xi2,yi2,zi2]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri2);

%% Computing slopes and slope direction
s=4;

xi1l=xi1(1:s:end,1:s:end);
yi1l=yi1(1:s:end,1:s:end);
zi1l=zi1(1:s:end,1:s:end);
xi2l=xi2(1:s:end,1:s:end);
yi2l=yi2(1:s:end,1:s:end);
zi2l=zi2(1:s:end,1:s:end);


[slope1,slope_dir1]=Slope(xi1l,yi1l,zi1l);
[slope2,slope_dir2]=Slope(xi2l,yi2l,zi2l);

figure; hold on;
surf(xi1(1:s:end,1:s:end),yi1(1:s:end,1:s:end),zi1(1:s:end,1:s:end));
StandardLight


figure; hold on;
surf(xi2(1:s:end,1:s:end),yi2(1:s:end,1:s:end),zi2(1:s:end,1:s:end));
StandardLight


AGUaxes
pcolorm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi,ri1(1:s:end,1:s:end));

AGUaxes
pcolorm(-fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi,ri2(1:s:end,1:s:end));

AGUaxes
pcolorm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi,slope1(1:s:end,1:s:end));

AGUaxes
pcolorm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi,slope2(1:s:end,1:s:end));

AGUaxes
pcolorm(fii1(1:s:end,1:s:end)/180*pi,lambdai1(1:s:end,1:s:end)/180*pi,...
    slope2(1:s:end,1:s:end)-slope1(1:s:end,1:s:end));

%% Plotting 

%
figure; hold on;
set(gca,'FontSize',12);

cond=(fii1(1:s:end,1:s:end)<50);

scatter((slope_dir1(cond)+slope_dir2(cond))/2*180/pi,slope1(cond)-slope2(cond),1,fii1(cond));
xlabel('Slope direction [deg]','FontSize',12);
ylabel('Slope difference [deg]','FontSize',12);
box on;

%
figure; hold on;
set(gca,'FontSize',12);
hist(slope1(:)-slope2(:),50);
xlabel('slope difference [deg]','FontSize',12);
ylabel('Number of points','FontSize',12);

%
% AGUaxes
% pcolorm(fii/180*pi,lambda/180*pi,slope1-slope2); shading interp;
% cbar=colorbar('FontSize',12);
% ylabel(cbar,'Slope difference','FontSize',12);






