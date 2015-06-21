ccc

lambda0=30/180*pi;
fi0=-90/180*pi;



fi=(-90:1:90)/180*pi;

lambda=(-180:1:180)/180*pi;

[fii,lambdai]=meshgrid(fi,lambda);


Cos1=sin(fii).*sin(fi0)+cos(fii).*cos(fi0).*cos(lambda0-lambdai);
Cos2=cos(fii-fi0).*cos(lambdai-lambda0);



AGUaxes; hold on
pcolorm(fii,lambdai,Cos1); shading interp;
[latc,lonc]=scircle1(fi0*180/pi,lambda0*180/pi,90);
plotm(latc/180*pi,lonc/180*pi,'--','LineWidth',5);
colorbar
plotm(fi0,lambda0,'*k','MarkerSize',40);
contourm(fii,lambdai,Cos1,-1:0.1:1,'-k','LineWidth',4);


AGUaxes; hold on
pcolorm(fii,lambdai,Cos2); shading interp;
plotm(latc/180*pi,lonc/180*pi,'--','LineWidth',5);
colorbar
plotm(fi0,lambda0,'*k','MarkerSize',40);
contourm(fii,lambdai,Cos2,-1:0.1:1,'-k','LineWidth',4);


AGUaxes; hold on
pcolorm(fii,lambdai,Cos1-Cos2); shading interp;
plotm(latc/180*pi,lonc/180*pi,'--','LineWidth',5);
colorbar
plotm(fi0,lambda0,'*k','MarkerSize',40);
contourm(fii,lambdai,Cos1-Cos2,-1:0.1:1,'-k','LineWidth',4);
