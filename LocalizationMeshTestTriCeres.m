ccc

%% figure settings

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Input Paramerers;

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150828_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150828_512.bds';

[~,shapename,~] = fileparts(shape_filename) ;

full_filename = [shape_folder shape_filename];

MaxDegreeTopo=100;
Resolution=0.5;
L=15;
MinConcentration=0.80;
NTess=3;
circle_rad=30;
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
        
        p=polyfit(log10(l(gd)),log10(sdl(gd,i)),2);
        
        aq=p(1);
        bq=p(2);
        
        xq(i) = mean(log10(l(gd)));
        
        curv(i) = 2*aq/((1+(bq+2*aq*xq(i)).^2).^1.5);
        bta(i) = p(1);
              
%         p = polyfit(log10(l(gd)),log10(sdl_mean(gd,j)),2);
    end
    
    progressbar(j/numel(fii));
    bta_mean(j)=mean(bta);
    curv_mean(j)=mean(curv);
    
    if (NumberOfTapers>1)
        sdl_mean(:,j)=mean(sdl,2);
        
    else
        sdl_mean(:,j)=(sdl);
    end
end

%% Plot latitude vs spectrum curv
figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

scatter(fii',curv_mean,40,lambdai','filled');

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

scatter(lambdai,fii',40,curv_mean,'filled');

xlim([-180 180]);
ylim([-90 90]);

colorbar
%% Plot power spectum slope as map

Li=20;
lmcosi_b=xyz2plm(curv_mean,Li,'irr',fii,lambdai);

[bta_mean_sh,lon,lat]=plm2xyz(lmcosi_b,0.5);
[lon,lat]=meshgrid(lon,lat);

AGUaxes;
pcolorm(lat,lon,bta_mean_sh); shading interp;

AGUaxes;
pcolorm(fi_grid*180/pi,lambda_grid*180/pi,Hi); shading interp;
colormap jet;
caxis([-7 7]);
colorbar

%% Clustering 

Y = log10(sdl_mean(gd,:));

NClucters=3;

% Admittance clusters

[idx,ctrs,sumd,D] = kmeans(Y',NClucters,'distance','sqEuclidean','onlinephase','on');
cc=jet(NClucters);
Dm=log10(min(D'));
[madm,ixs]=sort(mean(ctrs,2));
[~,ixs2]=sort(ixs);
cc2=cc(ixs2,:);

for i=1:NClucters
    ctrs_std(:,i)=std(Y(:,idx==i),0,2);
end


%%

figure; hold on;
set(gca,'FontSize',20);

for i=1:NClucters    
    plot(gd,ctrs(i,:),'-o','MarkerSize',3,'LineWidth',2,'Color',cc2(i,:));
    h=errorbar(gd,ctrs(i,:),ctrs_std(:,i),'LineWidth',1,'Color',cc2(i,:));
end

xlabel('SH Degree','FontSize',20);
ylabel('Topography Power []','FontSize',20);
box on;
grid on;

%% 

AGUaxes;
scatterm(fii,lambdai,100,idx,'filled');

% Voronoi diagram

xu=cos(fii/180*pi).*cos(lambdai/180*pi);
yu=cos(fii/180*pi).*sin(lambdai/180*pi);
zu=sin(fii/180*pi);

xyzu=[xu yu zu]';

[P, K, voronoiboundary] = voronoisphere(xyzu,'resolution', 0.1/180*pi);

% Graphic
% 
% f = figure(1);
% clf(f);
% set(f,'Renderer','zbuffer');
% ax = axes('Parent', f);
% hold(ax, 'on');
% axis(ax,'equal');
% 
% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');
% clmap = cool();
% ncl = size(clmap,1);
% 
% for k = 1:numel(xu)
%     X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
%     fill3(X(1,:),X(2,:),X(3,:),cc(idx(k),:),'Parent',ax,'EdgeColor','w');
% end
% 
% alpha(0.5);
% 
% axis(ax,'equal');
% axis(ax,[-1 1 -1 1 -1 1]);

%% Mapping hexagons

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    plotm(fi_b*180/pi,lambda_b*180/pi,'-r','LineWidth',2);
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end

AGUaxes; hold on
% plot3(ax, xyzu(1,:),xyzu(2,:),xyzu(3,:),'wo');
clmap = cool();
ncl = size(clmap,1);

for k = 1:numel(xu)
    X = voronoiboundary{k};
%     cl = clmap(mod(k,ncl)+1,:);
    [lambda_b,fi_b,~]=cart2sph(X(1,:),X(2,:),X(3,:));
    fillm(fi_b*180/pi,lambda_b*180/pi,0,cc2(idx(k),:),'EdgeColor','none');
%      fillm(fi_b,lambda_b,0,cc3(idc(k),:),'EdgeColor','none');
end

