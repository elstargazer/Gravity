ccc

%% figure settings

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Input Paramerers;

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPG/Survey_20151016/';
shape_filename='global.bds';

[~,shapename,~] = fileparts(shape_filename) ;
full_filename = [shape_folder shape_filename];

MaxDegreeTopo    = 100;
Resolution       = 1.0;
L                = 20;
MinConcentration = 0.85;
NTess            = 3;
circle_rad       = 20;

gd=4+L:MaxDegreeTopo-L+1;

a=481.000;
c=446.000;

%% Load topography model

[x_grid,y_grid,z_grid]=ReadSPC(full_filename,Resolution,'grid');

ri=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);
ri=flipud(ri');

% expand in SH
[~,~,Hi]=XYZ2BLH(x_grid,y_grid,z_grid,a,Eccentricity(a,c));
[lambda_grid,fi_grid,~]=cart2sph(x_grid,y_grid,z_grid);
lmcosi_shape= xyz2plm(ri,MaxDegreeTopo);

% check SH expansion by expanding back to spacial domain
[ri_sh,lon_sh,lat_sh] = plm2xyz(lmcosi_shape);
[lon_sh,lat_sh] = meshgrid(lon_sh,lat_sh);

AGUaxes;
pcolorm(lat_sh,lon_sh,ri_sh);

AGUaxes;
pcolorm(fi_grid*180/pi,lambda_grid*180/pi,Hi);

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

%% Write small circle coordinates

for i=1:numel(fii)
    
    [lat_sc,lon_sc] = scircle1(fii(i),lambdai(i),circle_rad);
    filename = [num2str(i) '.cl'];
    in = fopen(filename,'w');
    fprintf(in,'%6.2f %6.2f\n',[lon_sc lat_sc]');
    fclose(in);
    
end

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
        bta(i) = p(1);
        
        % p = polyfit(log10(l(gd)),log10(sdl_mean(gd,j)),2);
    end
    
    progressbar(j/numel(fii));
    
    if (NumberOfTapers>1)
        sdl_mean(:,j)=mean(sdl,2);
        sdl_std(:,j) = std(sdl,0,2);
        
    else
        sdl_mean(:,j)=(sdl);
        sdl_std(:,j) = NaN;
        
    end
end

%% Plot all power spectra

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'YScale','log');
xlim([L+3 MaxDegreeTopo-L]);

xlabel('Latitude [deg]','FontSize',fntsize);
ylabel('Spectral density','FontSize',fntsize);

plot(repmat(l,[1 size(sdl_mean,2)]),sdl_mean,'-','MarkerSize',2,'Color','k');

%% Plot power spectum slope as map

Li=23;
L_expand = 8;
sdl_pl = sdl_mean(find(l==Li),:);
lmcosi_pl=xyz2plm(sdl_pl,L_expand,'irr',fii,lambdai);

[pl,lon,lat]=plm2xyz(lmcosi_pl,0.5);
[lon,lat]=meshgrid(lon,lat);

AGUaxes;
pcolorm(lat,lon,pl); shading interp;

%% Power spectrum as a func of latitude

Li = L+3:1:MaxDegreeTopo-L;

i = 1;

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
xlim([-90 90]);

set(gca,'XTick',-90:30:90);

fi_lin = -90:1:90;

xlabel('Latitude [deg]','FontSize',fntsize);
h_ylab = ylabel(['Spectral density [km^2] at n = ' num2str(Li(1))],'FontSize',fntsize);

sdl_pl = sdl_mean(Li(i),:);
sdl_err_pl = sdl_std(Li(i),:);

pl_pow = plot(fii,sdl_pl,'.','MarkerSize',10,'Color','k');
% errorbar(fii,sdl_pl,sdl_err_pl,'.','MarkerSize',10,'Color','r');

[pfit,S] = polyfit(fii,sdl_pl',4);
[sdl_lin,sdl_lin_std] = polyconf(pfit,fi_lin,S);

% p1 = plot(fi_lin,sdl_lin,'-k','LineWidth',1);
% p2 = plot(fi_lin,sdl_lin-sdl_lin_std,'--r','LineWidth',1);
% p3 = plot(fi_lin,sdl_lin+sdl_lin_std,'--r','LineWidth',1);

[binCenters,binMean,binStandardDev]=rebin(fii,sdl_pl',3);
errorbar(binCenters,binMean,binStandardDev,'.r','LineWidth',2);

set(gca,'YLim',[max(0,min(sdl_lin-sdl_lin_std))...
    max(sdl_lin+sdl_lin_std)]);

which_deg_to_print = 23;

% waitforbuttonpress;

if Li(i) == which_deg_to_print
    disp(['printing for n = ' num2str(Li(i))]);
    PrintWhite([fig_folder 'Fig_LocalizedTopo_' ...
        num2str(Li(i)) '.jpg']);
end

which_deg_to_print = 40;


for i=2:numel(Li)
    sdl_pl = sdl_mean(Li(i),:);
    set(pl_pow,'YData',sdl_pl);
    h_ylab = ylabel(['Topography power at n = ' num2str(Li(i))],'FontSize',fntsize);
    
    [pfit,S] = polyfit(fii,sdl_pl',4);
    [sdl_lin,sdl_lin_std] = polyconf(pfit,fi_lin,S);
    
    %     set(p1,'YData',sdl_lin);
    %     set(p2,'YData',sdl_lin-sdl_lin_std);
    %     set(p3,'YData',sdl_lin+sdl_lin_std);
    
    set(gca,'YLim',[max(0,min(sdl_lin-sdl_lin_std))...
        max(sdl_lin+sdl_lin_std)]);
    
    
    [binCenters,binMean,binStandardDev]=rebin(fii,sdl_pl',3);
    errorbar(binCenters,binMean,binStandardDev,'.r','LineWidth',2);
       
    if Li(i) == which_deg_to_print
        disp(['printing for n = ' num2str(Li(i))]);
        PrintWhite([fig_folder 'Fig_LocalizedTopo_' ...
            num2str(Li(i)) '.jpg']);
    end
    
    waitforbuttonpress;
end

%%

figure('Color','k');
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize,'Color','k',...
    'YColor','w','XColor','w');
set(gca,'FontSize',fntsize);

hold on;grid on; box on;

i1 = find(Li==23);

sdl_pl = sdl_mean(Li(i1),:);
sdl_err_pl = sdl_std(Li(i1),:);

plot(fii,sdl_pl,'.','MarkerSize',10,'Color','w');
[binCenters,binMean,binStandardDev]=rebin(fii,sdl_pl',3);
errorbar(binCenters,binMean,binStandardDev,'.r','LineWidth',2);

xlabel('Latitude [deg]','FontSize',fntsize,'Color','w');
ylabel(['Spectral density [km^2] at n = ' num2str(Li(i1))],'FontSize',fntsize,'Color','w');

xlim([-90 90]);
set(gca,'XTick',-90:30:90);


%% Figure with subplots
figure;
set(gcf, 'Units','centimeters', 'Position',[0 0 26 9])
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
subplot(1,2,1);
hold on;grid on; box on;

i1 = find(Li==23);
i2 = find(Li==40);

sdl_pl = sdl_mean(Li(i1),:);
sdl_err_pl = sdl_std(Li(i1),:);

plot(fii,sdl_pl,'.','MarkerSize',10,'Color','k');
[binCenters,binMean,binStandardDev]=rebin(fii,sdl_pl',3);
errorbar(binCenters,binMean,binStandardDev,'.r','LineWidth',2);

xlabel('Latitude [deg]','FontSize',fntsize);
ylabel(['Spectral density [km^2] at n = ' num2str(Li(i1))],'FontSize',fntsize);

xlim([-90 90]);
set(gca,'XTick',-90:30:90);

subplot(1,2,2);
hold on;grid on; box on;

sdl_pl = sdl_mean(Li(i2),:);
sdl_err_pl = sdl_std(Li(i2),:);

plot(fii,sdl_pl,'.','MarkerSize',10,'Color','k');
[binCenters,binMean,binStandardDev]=rebin(fii,sdl_pl',3);
errorbar(binCenters,binMean,binStandardDev,'.r','LineWidth',2);

xlabel('Latitude [deg]','FontSize',fntsize);
ylabel(['Spectral density [km^2] at n = ' num2str(Li(i2))],'FontSize',fntsize);

xlim([-90 90]);
set(gca,'XTick',-90:30:90);

PrintWhite([fig_folder 'Fig_LocalizedTopo_23_40.jpg'])


%% Clustering

Y = log10(sdl_mean(gd,:));
NClucters=3;

%  spectrum clusters

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

%%
%
% figure;
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on;grid on; box on;
%
% set(gca,'YScale','log');
%
% step = 20;
% lat_center = -90:step:90;
%
% ccj = jet(numel(lat_center));
%
% for i=1:numel(lat_center)
%
%
%     ind = find(abs(fii - lat_center(i)) < step/2);
%
%     sdl_lat_avg = mean(sdl_mean(gd,ind),2);
%     sdl_lat_std =  std(sdl_mean(gd,ind),0,2);
%
%     errorbar(gd,sdl_lat_avg,sdl_lat_std,'Color',ccj(i,:));
%
% end
%
% xlabel('Degree','FontSize',fntsize);


%%

fi_lin = linspace(0,90,5);
ccj = jet(numel(fi_lin)-1);

ccj = [1 0.2 0.2;
       0.2 1 0.2;
       0.2 0.2 1;
       0 1 1];

figure('Color','k');
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize,'Color','k',...
    'YColor','w','XColor','w');
set(gca,'FontSize',fntsize);
set(gca,'YScale','Log');
set(gca,'XScale','Log');

hold on;box on;grid on;

leg_text = cell(1,numel(fi_lin)-1);
% leg_text{1} = 'global';
% plot(repmat(l,[1 numel(fii)]),sdl_mean,'-b');
% plot(l,sdl,'-k');

for i=1:numel(fi_lin)-1
    
    cond = (abs(fii) > fi_lin(i)) &  (abs(fii) < fi_lin(i+1));
    
    sdl_sub = sdl_mean(:,cond); 
    sdl_sub_mean = mean(sdl_sub,2);  
    sdl_sub_mean_std = std(sdl_sub,0,2);  
    
    n_spec = sum(cond);

    h(i) = plot(l(gd),sdl_sub_mean(gd),'-','Color',ccj(i,:),'LineWidth',3);
%     

if i==1
% h_err = errorbar(l(gd),sdl_sub_mean(gd),sdl_sub_mean_std(gd)/sqrt(n_spec),...
%          'Color',ccj(i,:),'LineWidth',1); 
     
     plot(l(gd),sdl_sub_mean(gd)+sdl_sub_mean_std(gd)/sqrt(n_spec),'.','Color',ccj(i,:),'MarkerSize',10);
     plot(l(gd),sdl_sub_mean(gd)-sdl_sub_mean_std(gd)/sqrt(n_spec),'.','Color',ccj(i,:),'MarkerSize',10);

end
    
      leg_text{i} = [num2str(fi_lin(i)) '^{\circ}' '- ' num2str(fi_lin(i+1)) '^{\circ}'];

end


legend(h,leg_text,'Color','k','TextColor','w','Color','k','TextColor','w');

% lmcosi_shape_m = lmcosi_shape;
% lmcosi_shape_m(:,3:4) = lmcosi_shape(:,3:4)*1000; 
[sdl,l] = plm2spec(lmcosi_shape);

xlabel('Spherical harmonic degree','FontSize',fntsize,'Color','w');
ylabel('Spectral density [km^{2}]','FontSize',fntsize,'Color','w');

PrintBlack([fig_folder 'Fig_LocalizedLatSpectra_k.jpg']);


















