ccc

%% figure settings

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Input Paramerers;

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150702_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150702_256.bds';

[~,shapename,~] = fileparts(shape_filename) ;

full_filename = [shape_folder shape_filename];

MaxDegreeTopo=35;
Resolution=0.75;
L=6;
MinConcentration=0.93;
NTess=3;
circle_rad=35;
gd=2+2+L:MaxDegreeTopo-L+1;

a=481.000;
c=446.000;

[x_grid,y_grid,z_grid]=ReadSPC(full_filename,Resolution,'grid');

ri=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);
ri=ri';

[~,~,Hi]=XYZ2BLH(x_grid,y_grid,z_grid,a,Eccentricity(a,c));
[lambda_grid,fi_grid,~]=cart2sph(x_grid,y_grid,z_grid);
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
        
        [sdl(:,i),l,~,~,~]=...
            plm2spec(lmcosi_shape_w{i});
        
        p=polyfit(log10(l(gd)),log10(sdl(gd,i)),1);
        bta(i)=p(1);
    end
    
    progressbar(j/numel(fii));
    bta_mean(j)=mean(bta);
    
    if (NumberOfTapers>1)
        sdl_mean(:,j)=mean(sdl,2);
    else
        sdl_mean(:,j)=(sdl);
    end
end

%% Plot latitude vs power spectrum slope
figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

scatter(fii',bta_mean,40,lambdai','filled');

xlabel('Latitude [deg]','FontSize',fntsize);
ylabel('Spectral slope','FontSize',fntsize);

cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Longitude [deg]','FontSize',fntsize);

xlim([-90 90]);
set(gca,'XTick',-90:30:90);

PrintWhite([fig_folder 'Fig_SpecSlope_' shapename '.jpg']);

%% Plot power spectum slope as scatter map

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

scatter(lambdai,fii',40,bta_mean,'filled');

xlim([-180 180]);
ylim([-90 90]);


%% Plot power spectum slope as map

Li=20;
lmcosi_b=xyz2plm(bta_mean,Li,'irr',fii,lambdai);

[bta_mean_sh,lon,lat]=plm2xyz(lmcosi_b,0.5);
[lon,lat]=meshgrid(lon,lat);

AGUaxes;
pcolorm(lat,lon,bta_mean_sh); shading interp;

AGUaxes;
pcolorm(fi_grid*180/pi,lambda_grid*180/pi,Hi); shading interp;
colormap jet;
caxis([-7 7]);
colorbar











