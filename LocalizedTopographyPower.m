ccc

%% Load topography model

MaxDegreeTopo=300;

load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
lmcosi_shape(4,3)=0;
[ri,lon,lat,~]=plm2xyz(lmcosi_shape,1,[],[],[]);

[loni,lati]=meshgrid(lon,lat);

MeanRadius=lmcosi_shape(1,3);

ri=ri-MeanRadius;

%% Geo grid

fi=-80:20:68;
lambda=0:20:340;

circle_rad=50;

[fii,lambdai]=meshgrid(fi,lambda);

L=5;
L2=100;

%% 
figure('Position',[1 1 1500 1500]);


sp2=subplot(2,1,2);

set(gca,'FontSize',20);

sp1=subplot(2,1,1);

ax=axesm('mollweid','frame','on','FontSize',24,'Grid','on','MLabelParallel',...
    'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
    , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',2,'FontColor',[0 0 0]...
    ,'FontSize',24,'GAltitude',Inf,'Geoid',[1 0]);

pcolorm(lati/180*pi,loni/180*pi,ri);


[latc,lonc] = scircle1(fii(1),lambdai(1),circle_rad);
pc=plotm(latc/180*pi,lonc/180*pi,'-b','LineWidth',5);



%% video

% Prepare the new file.
vidObj = VideoWriter('local_topo_power2.avi');
open(vidObj);



for j=1:numel(fii)

% patch coordinates
lat_center=fii(j);
lon_center=lambdai(j);

colat_center=90-lat_center;

% small circle

[latc,lonc] = scircle1(lat_center+0.1,lon_center+0.1,circle_rad);

delete(pc)

subplot(sp1)

pc=plotm(latc/180*pi,lonc/180*pi,'-b','LineWidth',5);



%% Using glmalphapto

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,lon_center,90-lat_center);

s=size(G2)


for i=1:s(2)
       
    lmcosi_window{i}=glm2lmcosi(G2,i);
    [w(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
end

%% Using localization

for i=1:s(2)
%    [r(:,:,i),lor,lar,Plm]=plm2xyz([jk1 jk2 C{i}],1);    
    
    

    [lmcosi_shape_w{i},~]=xyz2plm(w(:,:,i).*ri,L2,'im');
    
    Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});

   
end

% Z2_local=Sgt_local./Stt_local;

%% Localized admittance plot

Tapers=find(V2>0.7);

N0=(L+1)*(circle_rad/180);

% K=2;
% (K+1)*pi/(circle_rad/180*pi)-1

% plot_local=plot(Z_local(:,1:MaxTaper),'-.r');
% plot_local_shape=plot(Z_local_shape(:,1:MaxTaper),'--.r');

Stt_local_mean=mean(Stt_local(:,Tapers),2);

Stt_local_std=std(Stt_local(:,Tapers),0,2);

Degree=(1:numel(Stt_local_std))';


j/numel(fii)

subplot(sp2)

semilogy(Degree,Stt_local_mean,'LineWidth',2);

grid on;
ylabel('Topo power ','FontSize',20);
xlabel('Degree','FontSize',20);

% drawnow;
% pause(0.1)

currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);


end

% Close the file.
close(vidObj)