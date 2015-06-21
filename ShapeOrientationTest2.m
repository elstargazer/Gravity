ccc

%% Read Shape models
% ShapeModelFileName1='~/Dawn/Balmino/VestaTest/VestaPOSTVESTAshape_geoc_elev12m.grd';
% ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR12m.grd';

ShapeModelFileName1='~/Dawn/Balmino/VestaTest/Vesta20140513shape_geoc_elev_3m.grd';
ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR_3m.grd';

[lambdai1,fii1,ri1]=ReadGRD(ShapeModelFileName1); ri1=ri1*1000;
[lambdai2,fii2,ri2]=ReadGRD(ShapeModelFileName2);

[xi1,yi1,zi1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);
[xi2,yi2,zi2]=sph2cart(lambdai2/180*pi,fii2/180*pi,ri2);

a=280900;
c=226200;

[~,~,H1]=XYZ2BLH(xi1,yi1,zi1,a,Eccentricity(a,c));
[~,~,H2]=XYZ2BLH(xi2,yi2,zi2,a,Eccentricity(a,c));

clear xi1 yi1 zi1 xi2 yi2 zi2 ri2;

R=265000;

% Vesta20130522shape_geoc_elev_3m_gridline.grd
% VestaShapeDLR_3m_mean.grd 

%% Create centers

NTess=1;
r_circle=15;

TR=IcosahedronMesh;
TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;
x_t=TR_2.X(:,1);
y_t=TR_2.X(:,2);
z_t=TR_2.X(:,3);

[lambda_center,fi_center,~]=cart2sph(x_t,y_t,z_t);

lambda_center=lambda_center*180/pi;
fi_center=fi_center*180/pi;

%%
AGUaxes

for i=1:numel(fi_center)    
    [latc,lonc] = scircle1(fi_center(i),lambda_center(i),r_circle);
    plotm(latc/180*pi,lonc/180*pi,'-r','LineWidth',4);       
end

%% Cut 

ang_approx=[0 0 0];
angm=zeros(numel(lambda_center),3);

tic

progressbar(0);
% i=1;
parfor i=1:numel(fi_center)

    [dist,az] = distance(fii1,lambdai1,fi_center(i),lambda_center(i));
    Cond=(dist<r_circle);

    lambdac1=lambdai1(Cond)/180*pi;
    fic1=fii1(Cond)/180*pi;
    
    rc1=ri1(Cond);
%     rc2=ri2(Cond);

     Hc1=H1(Cond);
     Hc2=H2(Cond);
  
%     plotm(fic1,lambdac1,'.k','MarkerSize',1);

     options=optimset('TolX',10/3600/180*pi,'Display','iter');
    [angm(i,:),~] = fminsearch(@(ang) OrientationMetric(fic1,lambdac1,rc1,Hc1,Hc2,ang),ang_approx,options);
%     ang_approx=angm(i,:);
%     progressbar(i/numel(lambda_center));
     
end
progressbar(1);

t_el=toc

max_lat=50;

latcond=fi_center<max_lat;

angm(latcond,:)*180/pi*60
mean(angm(latcond,:)*180/pi*60)
std(angm(latcond,:)*180/pi*60)



