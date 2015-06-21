ccc

FileNameGEODYN='/Users/antonermakov/Dawn/Gravity/Frank20/grvfld.dawn.img.survlam5wim_5';
FileNameJPL='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

N_trunc=10;

%% Read MIT gravity model
[lmcosi,mu,Rref]=ReadGEODYNGravity(FileNameGEODYN);


% trancate 
lmcosi=TruncateGravityModel(lmcosi,N_trunc,1);


%% Read JPL gravity model
[lmcosi2,Rref2,mu2,~]=ReadGRAILGravityModel(FileNameJPL);


% truncate
lmcosi2=AddZeroHarm(lmcosi2,1);
lmcosi2=TruncateGravityModel(lmcosi2,N_trunc,1);


%% Plot errors

figure; hold on;

GravityModelSpectrum(lmcosi2,'r');
GravityModelSpectrum(lmcosi,'b')

legend({'JPL','JPL error','MIT','MIT error'})


%% Choose reference surface
arg=290000;
crg=265000;

AngleStep=1;

[fii,lambdai]=MakeGeoGrid(AngleStep,AngleStep);

xe=arg*cos(fii).*cos(lambdai);
ye=arg*cos(fii).*sin(lambdai);
ze=crg*sin(fii);

%% Compute gravity acceleration
[gx1,gy1,gz1]=GravityAcceleration(mu,Rref,lmcosi,xe,ye,ze);
[gx2,gy2,gz2]=GravityAcceleration(mu2,Rref2,lmcosi2,xe,ye,ze);

[gu1,~,~]=GravityComponents(gx1,gy1,gz1,xe,ye,ze,arg,crg);
[gu2,~,~]=GravityComponents(gx2,gy2,gz2,xe,ye,ze,arg,crg);

%% Map gravity acceleration
% 
MapRadialGrid(1e5*gu1');
title('grvfld.dawn.img.survlam5wim_5','FontSize',20);
MapRadialGrid(1e5*gu2');
title('JGV20G02.SHA','FontSize',20);

MapRadialGrid(1e5*(gu1-gu2)');
 
%% Correlation
figure;
set(gca,'FontSize',20);
Corr=SphericalHarmonicCorrelation(lmcosi,lmcosi2,'k','-o');
xlabel('Degree','FontSize',20);
ylabel('Correlation coefficient','FontSize',20);
xlim([2 20])
 
% [sdl1,l1,~,~,~,~]=plm2spec(lmcosi);
% [sdl2,l2,~,~,~,~]=plm2spec(lmcosi2);
% 
% lmcosid=lmcosi;
% lmcosid(:,3:4)=lmcosi(:,3:4)-lmcosi2(:,3:4);
% 
% [sdld,ld,~,~,~,~]=plm2spec(lmcosid);
% 
% figure
% hold on;
% plot(l1,log10(sdl1),'ro-','LineWidth',2);
% plot(l2,log10(sdl2),'go-','LineWidth',2);
% plot(ld,log10(sdld),'bo-','LineWidth',2);
% ylabel('log_{10}(Coeff Mag)','FontSize',20)
% xlabel('Degree','FontSize',20)
% 
% legend({'MIT5','JPL20','Diff'},'FontSize',20);
% 
