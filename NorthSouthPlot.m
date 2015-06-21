filename1='/Users/ermakov/Dawn/ShapeModel/Gaskell/20140513/Vesta20140513shape_geoc_elev.grd';
filename2='/Users/ermakov/Dawn/ShapeModel/DLR/DLR_SHAPE_2014_06_19/VestaShapeDLR64.grd';

[lambdai1,fii1,ri1]=ReadGRD(filename1);
[~       ,~   ,ri2]=ReadGRD(filename2); ri2=flipud(ri2/1000);

[xi1,yi1,zi1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);
[xi2,yi2,zi2]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri2);


% south
figure; hold on;
title('SPC South','FontSize',20);
surf(xi1(1:200,:),yi1(1:200,:),zi1(1:200,:));
lighting phong
shading interp
axis tight equal 
light('Position',[-1 1 -1],'Style','infinite');
xlabel('X [km]','FontSize',20);
ylabel('Y [km]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius [km]','FontSize',20);
view(270,-30)

figure; hold on;
title('SPG South','FontSize',20);
set(gca,'FontSize',20);
surf(xi2(1:200,:),yi2(1:200,:),zi2(1:200,:));
lighting phong
shading interp
axis tight equal 
light('Position',[-1 1 -1],'Style','infinite');
xlabel('X [km]','FontSize',20);
ylabel('Y [km]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius [km]','FontSize',20);
view(270,-30)

% north
figure; hold on;
title('SPC North','FontSize',20);
set(gca,'FontSize',20);
surf(xi1(end-200:end,:),yi1(end-200:end,:),zi1(end-200:end,:));
lighting phong
shading interp
axis tight equal 
light('Position',[-1 1 1],'Style','infinite');
xlabel('X [km]','FontSize',20);
ylabel('Y [km]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius [km]','FontSize',20);
view(270,30)

figure; hold on;
title('SPG North','FontSize',20);
set(gca,'FontSize',20);
surf(xi2(end-200:end,:),yi2(end-200:end,:),zi2(end-200:end,:));
lighting phong
shading interp
axis tight equal 
light('Position',[-1 1 1],'Style','infinite');
xlabel('X [km]','FontSize',20);
ylabel('Y [km]','FontSize',20);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius [km]','FontSize',20);
view(270,30)

