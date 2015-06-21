function ShowShapeWithGrid(lmcosi_shape,MaxDegreeTopo,Resolution,GridResolution,MaxLatitude)
 
HowMuchUp=1000;

%% Compute shape model
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);  
[ri,~]=plm2xyz(lmcosi_shape,Resolution);
s=size(ri);
fi=linspace(90,-90,s(1))/180*pi;
lambda=linspace(0,360,s(2))/180*pi;

[lambdai,fii]=meshgrid(lambda,fi);
[xi,yi,zi]=sph2cart(lambdai,fii,ri);
[xig,yig,zig]=sph2cart(lambdai,fii,ri+HowMuchUp);

%% Plot shape model
figure; hold on;
hs=surf(xi,yi,zi,ri); shading interp;
StandardLight
set(gca,'FontSize',20);
axis equal

%% Plot grid
fig=(90:-GridResolution:-90)/180*pi;
lambdag=(0:GridResolution:360)/180*pi;

[~, afi, ~] = setxor(fi,fig);
[~, alambda, ~] = setxor(lambda,lambdag);
fi(afi) = NaN;
lambda(alambda) = NaN;
[lambdaig,fiig] = meshgrid(lambda,fi);

MaxLatitude=MaxLatitude/180*pi;

idx2 = (~isnan(fiig)|~isnan(lambdaig))&(fii>-MaxLatitude)&(fii<MaxLatitude);

xg=nan(size(xig));
yg=nan(size(yig));
zg=nan(size(zig));

xg(idx2)=xig(idx2);
yg(idx2)=yig(idx2);
zg(idx2)=zig(idx2);

plot3(xg,yg,zg,'-k')
plot3(xg',yg',zg','-k')
