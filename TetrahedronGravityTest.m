% ccc

rr=301000;

rho=3000;

fi    =[1 1 2]/180*pi;
lambda=[0 2 1]/180*pi;


r=[300000 300000 300000];

[x,y,z]=sph2cart(lambda,fi,r);


FV=[1 2 3];


fir=(-90:1:90)/180*pi;
lambdar=(0:1:360)/180*pi;
[fir,lambdar]=meshgrid(fir,lambdar);
[xr,yr,zr]=sph2cart(lambdar,fir,rr);



%% Plot configuration

figure; hold on;
set(gca,'FontSize',20);
trisurf(FV,x,y,z);
plot3(0,0,0,'.');
plot3(xr,yr,zr,'.r');

%% Compute gravity 1


U1=GravityPotentialTriDen93(x,y,z,FV,xr,yr,zr,rho);
[ax1,ay1,az1]=GravityAccelerationTriDen(x,y,z,FV,xr,yr,zr,rho);

AGUaxes
set(gca,'FontSize',20);
pcolorm(fir,lambdar,U1);
plotm(fi,lambda,'o','MarkerEdgeColor','k','MarkerFaceColor','w');
cbar=colorbar('FontSize',20);
shading interp;


%% Compute gravity 2

xmean=sum(x)/4;
ymean=sum(y)/4;
zmean=sum(z)/4;

xm=x-xmean;
ym=y-ymean;
zm=z-zmean;

xm=[xm -xmean];
ym=[ym -ymean];
zm=[zm -zmean];

FVm=[ 1 2 3;
      1 4 2;
      2 4 3;
      3 4 1 ];
    


figure; hold on;
set(gca,'FontSize',20);
trisurf(FVm,xm,ym,zm);
plot3(0,0,0,'.');
alpha(0.5);
axis equal


xrm=xr-xmean;
yrm=yr-ymean;
zrm=zr-zmean;

plot3(xrm,yrm,zrm,'.r');


U2=GravityPotentialTri93(xm,ym,zm,FVm,xrm,yrm,zrm,rho);
[U2b,ax2,ay2,az2]=GravityPotAccTri(xm,ym,zm,FVm,xrm,yrm,zrm,rho);


AGUaxes
set(gca,'FontSize',20);
pcolorm(fir,lambdar,(U1-U2b)./U1);
plotm(fi,lambda,'o','MarkerEdgeColor','k','MarkerFaceColor','w');
cbar=colorbar('FontSize',20);
shading interp;

AGUaxes
set(gca,'FontSize',20);
pcolorm(fir,lambdar,U2./U1);
plotm(fi,lambda,'o','MarkerEdgeColor','k','MarkerFaceColor','w');
cbar=colorbar('FontSize',20);
shading interp;









