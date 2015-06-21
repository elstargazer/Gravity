ccc

%% Input Paramerers;
filename='/Users/antonermakov/Dawn/CeresShapeModel/OpNav/OpNav5/ceres_opnav5_512.bds';

MaxDegreeTopo=35;
Resolution=0.75;
L=5;
MinConcentration=0.93;
NTess=2;
circle_rad=35; 
GoodDegrees=2+L:MaxDegreeTopo-L;

[x_grid,y_grid,z_grid]=LoadOpNavShape(filename,Resolution,'grid');
ri=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);
ri=ri';

%% Load topography model

lmcosi_shape= xyz2plm(ri,MaxDegreeTopo);

%% Icosahedron mesh;
TR=IcosahedronMesh;
TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;
x_t=TR_2.X(:,1);
y_t=TR_2.X(:,2);
z_t=TR_2.X(:,3);

[lambdai,fii,~]=cart2sph(x_t,y_t,z_t);
lambdai=lambdai*180/pi;
fii=fii*180/pi;

% fii=fii(1:10);
Npoints=numel(fii);
% MaxDegreeExp=fix(0.5*(-3+sqrt(1+8*Npoints)))

%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,0,0);

figure; hold on;
set(gca,'FontSize',20);

plot(-sort(-V2),'-ok','LineWidth',3);

xlabel('Taper number','FontSize',20);
ylabel('Concentration factor [] ','FontSize',20);
box on;
grid on;
xlim([1 numel(V2)]);

s=size(G2);

for i=1:s(2)       
    lmcosi_window_basic{i}=glm2lmcosi(G2,i);    
end

Tapers=find(V2>MinConcentration);
NumberOfTapers=numel(Tapers);
lmcosi_window_basic=lmcosi_window_basic(Tapers);
NTapers=sum(V2>MinConcentration);


progressbar(0);

for j=1:numel(fii)

% patch coordinates
    lat_center=fii(j);
    lon_center=lambdai(j);

    colat_center=90-lat_center;
    alp=0;
    beta=lat_center-90; 
    gam=-lon_center;
    
    for i=1:NumberOfTapers
        
        [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,beta,gam,'dlmb');    
        [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},Resolution);    
        [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,MaxDegreeTopo,'im'); 
   
        [sdl(:,i),l,bta(i),lfit(:,i),logy(:,i)]=...
        plm2spec(lmcosi_shape_w{i}); 

    end    
    
    progressbar(j/numel(fii));
    bta_mean(j)=mean(bta);
end
progressbar(1);

figure; hold on;
set(gca,'FontSize',20);

plot(fii',bta_mean,'.');

xlabel('Latitude [deg]','FontSize',20);
ylabel('Spectral slope []','FontSize',20);

xlim([-90 90]);

figure; hold on;
set(gca,'FontSize',20);
scatter(lambdai',fii',100,bta_mean,'filled');
xlim([-180 180]);
ylim([-90 90]);

Li=10;
lmcosi_b=xyz2plm(bta_mean,Li,'irr',fii,lambdai);

[bta_mean_sh,lon,lat]=plm2xyz(lmcosi_b,0.5);

[lon,lat]=meshgrid(lon,lat);

AGUaxes;
pcolorm(lat,lon,bta_mean_sh); shading interp;














