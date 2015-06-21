ccc

% ShapeModelFileName1='~/Dawn/Balmino/VestaTest/VestaPOSTVESTAshape_geoc_elev12m.grd';
% ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR12m.grd';

% 
%ShapeModelFileName1='~/Dawn/Balmino/VestaTest/Vesta20130522shape_geoc_elev.grd';
% ShapeModelFileName2='~/Dawn/Balmino/VestaTest/VestaShapeDLR64.grd';

ShapeModelFileName1='~/Dawn/MATLAB/VestaShapeDLR64.grd';

ShapeModelList='ShapeModelList';

in=fopen(ShapeModelList,'r');

drnew_global=[];

[lambdai1,fii1,ri1]=ReadGRD(ShapeModelFileName1); 
% ri1=ri1*1000;
ri1=flipud(ri1);

fii1=fii1(1:4:end,1:4:end);
lambdai1=lambdai1(1:4:end,1:4:end);
ri1=ri1(1:4:end,1:4:end);

ShapeModelFileName2=fgetl(in);

j2off_center_global=[];
j2off_mean_global=[];
j2off_std_global=[];

while (ShapeModelFileName2~=-1)
    
ShapeModelFileName2='~/Dawn/MATLAB/Vesta20140513shape_geoc_elev.grd';
 
[lambdai2,fii2,ri2]=ReadGRD(ShapeModelFileName2);
ri2=ri2*1000;

fii2=fii2(1:4:end,1:4:end);
lambdai2=lambdai2(1:4:end,1:4:end);
ri2=ri2(1:4:end,1:4:end);


% ri2=flipud(ri2);
% ri2=circshift(ri2,[0 8640]);

resm=fii2(2,1)-fii2(1,1);
nf=size(ri2,1);
nl=size(ri2,2);


skip=10;

% AGUaxes; set(gca,'FontSize',20);
% pcolorm(fii1(1:skip:end,1:skip:end)/180*pi,...
%     lambdai1(1:skip:end,1:skip:end)/180*pi,....
%     ri1(1:skip:end,1:skip:end));
% colorbar('FontSize',20);
% 
% AGUaxes; set(gca,'FontSize',20);
% pcolorm(fii2(1:skip:end,1:skip:end)/180*pi,...
%     lambdai2(1:skip:end,1:skip:end)/180*pi,...
%     ri2(1:skip:end,1:skip:end));
% colorbar('FontSize',20);


% ri2=flipud(ri2);

%% Plotting difference

dr=ri2-ri1;
% 
% AGUaxes; set(gca,'FontSize',20);
% 
% pcolorm(fii2(1:skip:end,1:skip:end)/180*pi,...
%     lambdai2(1:skip:end,1:skip:end)/180*pi,...
%           dr(1:skip:end,1:skip:end))
%       
%   
% caxis([-400 400])
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation difference [m] ','FontSize',20);


% f1=figure; hold on;
% set(gca,'FontSize',20);
% hist(log10(abs(dr(:))),50,'FaceColor','b','EdgeColor','w')
% xlabel('log_{10}(\Delta r)','FontSize',20);
% ylabel('# of points','FontSize',20);
% box on;

% figure; hold on;
% plot(ri2(450,:))
% plot(ri1(450,:),'r')



% %% Lon. aver. rad
% rlonav=mean(ri1,2);
% rstd=std(ri1,0,2);
% 
% figure; hold on;
% set(gca,'FontSize',20);
% plot(fii1(:,1),rlonav,'r.','MarkerSize',1)
% 
% % plot(binEdges,binMean*60)
% % h=errorbar(binEdges,binMean*60,binStandardDev*60,'LineWidth',2);
% % errorbar_tick(h,1000);
% 
% xlabel('Latitude [deg]','FontSize',20);
% ylabel('Longitude averaged radius','FontSize',20);
% set(gca,'XTick',-90:30:90);
% xlim([-90 90]);
% % ylim([0 2000e6]);
% box on;


%%
row=450;
%WriteXYZ(lambdai1,fii1,dr,'DLRmGaskell12min.xyz');


% [x1,y1,z1]=sph2cart(lambdai1(1:10:end,1:10:end)/180*pi,fii1(1:10:end,1:10:end)/180*pi,ri1(1:10:end,1:10:end));
% [x2,y2,z2]=sph2cart(lambdai2(1:10:end,1:10:end)/180*pi,fii2(1:10:end,1:10:end)/180*pi,ri2(1:10:end,1:10:end));
% 
% figure
% surf(x1,y1,z1,ri1(1:10:end,1:10:end)); shading interp; StandardLight; axis equal;
% % caxis([210000 293000])
% 
% figure
% surf(x2,y2,z2,ri2(1:10:end,1:10:end)); shading interp; StandardLight; axis equal;
% caxis([210000 293000])


%% Searching Z-angle offset
exclude=10;
min_angle=-5;
max_angle=5;
step_angle=0.05;

[ri1new,delta,angle_min]=FindLongitudeOffset(lambdai1,ri1,ri2,exclude,...
    min_angle,max_angle,step_angle);

[~,delta_check,angle_min_check]=FindLongitudeOffset(lambdai1,ri1new,ri2,exclude,...
    min_angle,max_angle,step_angle);


[binEdges,binMean,binStandardDev]=rebin(fii1(exclude:(size(ri1,1)-exclude),1),angle_min(:),3);
% 
figure('Color','k'); hold on;
set(gca,'FontSize',20,'Color','k');
plot(fii1(exclude:(size(ri1,1)-exclude),1),-angle_min*60,'r.','MarkerSize',1)

plot(binEdges,-binMean*60)
h=errorbar(binEdges,-binMean*60,binStandardDev*60,'LineWidth',2);
errorbar_tick(h,1000);

xlabel('Latitude [deg]','FontSize',20,'Color','w');
ylabel('Longutude offset [min]','FontSize',20,'Color','w');
set(gca,'XTick',-90:30:90);
xlim([-90 90]);
ylim([min_angle max_angle]);

set(gca,'YColor','w') 
set(gca,'XColor','w') 

box on;

ylim([-2 2]*4);



%% Modifying grid

drnew=ri2-ri1new;

% AGUaxes; set(gca,'FontSize',20);
% pcolorm(fii2(1:skip:end,1:skip:end)/180*pi,...
%     lambdai2(1:skip:end,1:skip:end)/180*pi,...
%     drnew(1:skip:end,1:skip:end));
% caxis([-400 400])
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation difference [m] ','FontSize',20);
% WriteXYZ(lambdai1,fii1,drnew,'DLRmGaskell12min_mod.xyz');

% figure(f1)
% hist(log10(abs(drnew(:))),50,'FaceColor','r','EdgeColor','w')

%% j2 difference

[ri1new2]=ZonalDifference(ri1,ri2);


[binEdges,binMean,binStandardDev]=rebin(fii1(:),drnew(:),1);

figure('Color','k'); hold on;
set(gca,'FontSize',20);

set(gca,'YColor','w') 
set(gca,'XColor','w') 
set(gca,'Color','k');

box on;

% errorbar(fii1(:,1),drnewlonav,drnewstd);

plot(binEdges,-binMean,'k','LineWidth',3)
plot(fii1(:,1),-drnew(:,1:10:end),'.r','MarkerSize',1)
h=errorbar(binEdges,-binMean,binStandardDev,'LineWidth',1);
errorbar_tick(h,1000);

xlabel('Latitude [deg]','FontSize',20,'Color','w');
ylabel('Longitude-averaved difference [m] ','FontSize',20,'Color','w');
xlim([-90 90]);
ylim([-50 250]);
set(gca,'XTick',-90:30:90);
box on;
% 
% haxes1=gca;
% haxes1_pos = get(haxes1,'Position'); % store position of first axes
% 
% haxes2 = axes('Position',haxes1_pos,...
%               'XAxisLocation','top',...
%               'YAxisLocation','right',...
%               'Color','none');
%           
%           
% plot(haxes2,fii1(:,1),rlonav,'-k','LineWidth',3)
% line(fii1(:,1),rlonav/1000,'Parent',haxes2,'Color','k','LineWidth',4)
% set(haxes2,'xlim',[-90 90])
% set(haxes2,'ylim',[210 305])
% set(haxes2,'FontSize',20);
% ylabel(haxes2,'Radius [km]','FontSize',20);
% set(haxes2,'xtick',-90:30:90);

% 
% AGUaxes; set(gca,'FontSize',20);
% pcolorm(fii2(1:skip:end,1:skip:end)/180*pi,...
%     lambdai2(1:skip:end,1:skip:end)/180*pi,...
%     drnew2(1:skip:end,1:skip:end));
% caxis([-400 400])
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Elevation difference [m] ','FontSize',20);

drnew2=ri2-ri1new2;

j2off_center_global=[j2off_center_global; binEdges];
j2off_mean_global=[j2off_mean_global; binMean];
j2off_std_global=[j2off_std_global; binStandardDev];

ShapeModelFileName2=fgetl(in)

end



figure; hold on;
set(gca,'FontSize',20);
% errorbar(fii1(:,1),drnewlonav,drnewstd);

cc2=jet(size(j2off_center_global,1));


xlabel('Latitude [deg]','FontSize',20);
ylabel('Longitude-averaved difference [m] ','FontSize',20);
xlim([-90 90]);
ylim([-250 250]);
set(gca,'XTick',-90:30:90);
box on;


for i=1:size(j2off_center_global,1)

plot(j2off_center_global(i,:),j2off_mean_global(i,:),'LineWidth',3,'Color',cc2(i,:))
h=errorbar(j2off_center_global(i,:),j2off_mean_global(i,:),j2off_std_global(i,:),...
    'LineWidth',1,'Color',cc2(i,:));
errorbar_tick(h,1000);

end


% WriteXYZ(lambdai1,fii1,drnew,'DLRmGaskell12min_mod2.xyz');
% WriteXYZ(lambdai2,fii2,drnew2,'SPGminusSPC_3m.xyz');

% figure(f1)
% hist(log10(abs(drnew2(:))),50,'FaceColor','g','EdgeColor','w');


% %% SH
% MaxDeg=200;

% lmcosi_g=xyz2plm(ri1,MaxDeg);
% lmcosi_gnew2=xyz2plm(ri1new2,MaxDeg);
% lmcosi_d=xyz2plm(ri2,MaxDeg);

% Corr=SphericalHarmonicCorrelation(lmcosi_d,lmcosi_g);
% Corrnew2=SphericalHarmonicCorrelation(lmcosi_d,lmcosi_gnew2);
% 
% figure('Color','k'); hold on;
% set(gca,'FontSize',20,'Color','k');
% 
% ylabel('Correlation []','FontSize',20);
% 
% plot(Corr,'k-','LineWidth',1);
% plot(Corrnew2,'r-','LineWidth',1);

%% Histograms


[dr_db,~]=DebiasHist(dr,fii1);
[drnew_db,~]=DebiasHist(drnew,fii1);
[drnew2_db,fii1_db]=DebiasHist(drnew2,fii1);




% log1=log10(abs(dr(:)));
% log2=log10(abs(drnew(:)));
% log3=log10(abs(drnew2(:)));

log1=log10(abs(dr_db(:)));
log2=log10(abs(drnew_db(:)));
log3=log10(abs(drnew2_db(:)));

log1(log1<-20)=-2;
log2(log2<-20)=-2;
log3(log3<-20)=-2;

lat=fii1_db(:);


figure('Color','k'); hold on;
set(gca,'FontSize',20,'Color','k');
grid on;

set(gca,'YColor','w') 
set(gca,'XColor','w') 

xlabel('\Delta r [m] ','FontSize',20,'Color','w');
ylabel('# of points','FontSize',20,'Color','w');
box on;
xlim([-2 4]);
ylim([0 1.7e6]);


[n1, xout1] = hist(log1,50);
bar(xout1,n1,'r'); 

[n2, xout2] = hist(log2(~isnan(log2)),50);
bar(xout2,n2,'g');

[n3, xout3] = hist(log3(~isnan(log3)),50);
bar(xout3,n3,'b');

[n4, xout4] = hist(log3(lat<50),50);
bar(xout4,n4,'y');

alpha(0.7)

% legend({'Uncorrected '},'FontSize',20,'Color','w');
% legend({'Uncorrected ','Lon. offset corrected '},'FontSize',20);
% legend({'Uncorrected ','Lon. offset corrected ','J2 corrected '},'FontSize',20);
legend({'Uncorrected ','Lon. offset corrected ','J2 corrected ','Below 50 N   '},...
    'FontSize',20,'Color','w');

xtick=get(gca,'XTickLabels')
set(gca,'XTickLabels',{'0.01','0.1','1','10','100','1000','10000'})

%% Slopes

[xi1,yi1,zi1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1);
[xi1new,yi1new,zi1new]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1new);
[xi1new2,yi1new2,zi1new2]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1new2);

[xi2,yi2,zi2]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri2);

angle=zeros(1,numel(xi1));
anglenew=zeros(1,numel(xi1));
anglenew2=zeros(1,numel(xi1));

[nxi1,nyi1,nzi1] = surfnorm(xi1,yi1,zi1);
[nxi1new,nyi1new,nzi1new] = surfnorm(xi1new,yi1new,zi1new);
[nxi1new2,nyi1new2,nzi1new2] = surfnorm(xi1new2,yi1new2,zi1new2);

[nxi2,nyi2,nzi2] = surfnorm(xi2,yi2,zi2);

parfor i=1:numel(xi1)
    
    rvec1=[nxi1(i) nyi1(i) nzi1(i)];
    rvec1new=[nxi1new(i) nyi1new(i) nzi1new(i)];
    rvec1new2=[nxi1new2(i) nyi1new2(i) nzi1new2(i)];
    
    rvec2=[nxi2(i) nyi2(i) nzi2(i)];
    
    angle(i) = atan2(norm(cross(rvec1,rvec2)),dot(rvec1,rvec2));   
    anglenew(i) = atan2(norm(cross(rvec1new,rvec2)),dot(rvec1new,rvec2));   
    anglenew2(i) = atan2(norm(cross(rvec1new2,rvec2)),dot(rvec1new2,rvec2));   
end

angle=reshape(angle,size(xi1));
anglenew=reshape(anglenew,size(xi1));
anglenew2=reshape(anglenew2,size(xi1));

   
[angle_db,~]=DebiasHist(angle,fii1);
[anglenew_db,~]=DebiasHist(anglenew,fii1);
[anglenew2_db,fii1_db]=DebiasHist(anglenew2,fii1);
 
logs1=log10(abs(angle_db*180/pi));
logs2=log10(abs(anglenew_db*180/pi));
logs3=log10(abs(anglenew2_db*180/pi));

figure('Color','k'); hold on;
set(gca,'FontSize',20,'Color','k');
grid on;

logs1(logs1<-4)=-4;
logs2(logs2<-4)=-4;
logs3(logs3<-4)=-4;


[n1, xout1] = hist(logs1(~isnan(logs1)),50);
bar(xout1,n1,'r');

[n2, xout2] = hist(logs2(~isnan(logs2)),50);
bar(xout2,n2,'g');

[n3, xout3] = hist(logs3(~isnan(logs3)),50);
bar(xout3,n3,'b');

[n4, xout4] = hist(logs3(lat<50),50);
bar(xout4,n4,'y');


set(gca,'YColor','w') 
set(gca,'XColor','w') 

xlabel('\Delta(slope) [deg] ','FontSize',20,'Color','w');
ylabel('# of points','FontSize',20,'Color','w');
box on;
xlim([-1 1.5]);
% ylim([0 1.7e6]);

legend({'Uncorrected ','Lon. offset corrected ','J2 corrected ','Below 50 deg N   '},...
    'FontSize',20,'Color','w');

set(gca,'XTickLabel',{'0.1','0.3','1','3','10','30'})







%% COF offset

% dx=0;
% dy=0;
% dz=-0;
% 
% [x1,y1,z1]=sph2cart(lambdai1/180*pi,fii1/180*pi,ri1new2);
% 
% x1=x1+dx;
% y1=y1+dy;
% z1=z1+dz;
% 
% 
% [lambdai1n,fii1n,ri1new3]=cart2sph(x1,y1,z1);
% 
% 
% drnew3=ri2-ri1new3;
% 
% AGUaxes; set(gca,'FontSize',20);
% pcolorm(fii1/180*pi,lambdai1/180*pi,drnew3);
% caxis([-400 400])
% colorbar('FontSize',20);


% ri_g=plm2xyz(lmcosi_g,1);
% ri_d=plm2xyz(lmcosi_d,1);
% 
% Rref=265000;
% mu=17.28e9;
% 
% lmcosi_grav_d=TopoSH2GravitySH(ri1,mu,3456,Rref,200,20,5);
% lmcosi_grav_g=TopoSH2GravitySH(ri2,mu,3456,Rref,200,20,5);
% 
% zcm_d=lmcosi_grav_d(2,3)*Rref/NormCoef(1,0);
% xcm_d=lmcosi_grav_d(3,3)*Rref/NormCoef(1,1);
% ycm_d=lmcosi_grav_d(3,4)*Rref/NormCoef(1,1);
% 
% [xcm_d ycm_d zcm_d]
% 
% zcm_g=lmcosi_grav_g(2,3)*Rref/NormCoef(1,0);
% xcm_g=lmcosi_grav_g(3,3)*Rref/NormCoef(1,1);
% ycm_g=lmcosi_grav_g(3,4)*Rref/NormCoef(1,1);
% 
% [xcm_d ycm_d zcm_d]
% 
% [xcm_g ycm_g zcm_g]
% 
% [xcm_g ycm_g zcm_g]-[xcm_d ycm_d zcm_d]
% 
% 
% C20_g=lmcosi_grav_g(4,3);
% C20_d=lmcosi_grav_d(4,3);
% % 
% C21_g=lmcosi_grav_g(5,3);
% C21_d=lmcosi_grav_d(5,3);
% % 
% S21_g=lmcosi_grav_g(5,4);
% S21_d=lmcosi_grav_d(5,4);
% % 
% C22_g=lmcosi_grav_g(6,3);
% C22_d=lmcosi_grav_d(6,3);
% % 
% S22_g=lmcosi_grav_g(6,4);
% S22_d=lmcosi_grav_d(6,4);
% 
% 
% 
% lambda_eq_d=atan2(S22_d,C22_d)*180/pi;
% lambda_eq_g=atan2(S22_g,C22_g)*180/pi;
% 
% 
% (lambda_eq_d-lambda_eq_g)*60







