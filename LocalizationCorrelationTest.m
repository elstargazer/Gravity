%% Input Paramerers;
GravityFileName='/Users/antonermakov/GRAIL/Gravity_Models/jggrx_0660b_sha.tab.txt';
GravityFromTopoFileName='/Users/antonermakov/GRAIL/Gravity_Models/jggrx_0660b_sha.tab.txt';
% GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

MaxDegreeTopo=1800;
MaxDegreeGrav=500;
% Rref_gravity=293000;
Rref_gravity=1738000;
rhomean=3344;
Resolution=.1;
MaxTopoPower=4;
L=20;
MinConcentration=0.93;
NTess=2;
circle_rad=10; 
MaxDegreeExp=7;
rho_mean=3457.5;
GoodDegrees=2+L:MaxDegreeGrav-L;

%% Load gravity model

[lmcosi_grav,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

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

Z_t=zeros(numel(GoodDegrees),numel(fii));
Z_t_std=zeros(numel(GoodDegrees),numel(fii));
C_t=zeros(numel(GoodDegrees),numel(fii));
C_t_std=zeros(numel(GoodDegrees),numel(fii));
Z_tm=zeros(numel(GoodDegrees),NTapers*numel(fii));
C_tm=zeros(numel(GoodDegrees),NTapers*numel(fii));

progressbar(0);

for j=1:numel(fii)

% patch coordinates
    lat_center=fii(j);
    lon_center=lambdai(j);

    colat_center=90-lat_center;
    alp=0;
    bta=lat_center-90; 
    gam=-lon_center;
    
%     [fig_cm,ax_cm]=AGUaxes;
%     set(fig_cm,'Position',[1 1 1000*1.6 500*1.6]);
%     vidObj = VideoWriter('Tapers_Moon.avi');
%     vidObj.FrameRate=5;
%     vidObj.Quality=100;
%     open(vidObj);

    for i=1:NumberOfTapers
        
        [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,bta,gam,'dlmb');    
        [r(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},Resolution);    
    
        [lmcosi_fa_w{i},~]=xyz2plm(r(:,:,i).*fa,MaxDegreeGrav,'im');    
        [lmcosi_shape_w{i},~]=xyz2plm(r(:,:,i).*ri,MaxDegreeGrav,'im');
        
        
        %plotting tapers
        [lori,lari]=meshgrid(lor/180*pi,lar/180*pi);
        
        
%         AGUaxes
%         [latcirc,loncirc]=scircle1(lat_center,lon_center,circle_rad);
%         pcolorm(lari,lori,r(:,:,i)); shading interp;
%         plotm(latcirc/180*pi,loncirc/180*pi,'-k','LineWidth',3)
%         caxis([-1.5 1.5]);
%         title(['Taper # ' num2str(i) ' '],'FontSize',24);
        
%         lmcosi_shape_w{i}(:,3:4)=lmcosi_shape_w{i}(:,3:4)./lmcosi_shape_w{i}(1,3);
%         lmcosi_fa_w{i}(:,3:4)=lmcosi_fa_w{i}(:,3:4)./lmcosi_fa_w{i}(1,3);
         
        Sgt_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_shape_w{i});
        Stt_local(:,i)=CrossPower(lmcosi_shape_w{i},lmcosi_shape_w{i});    
        Sgg_local(:,i)=CrossPower(lmcosi_fa_w{i},lmcosi_fa_w{i});
            
%          currFrame = getframe(gcf);
%          writeVideo(vidObj,currFrame);

    end   
    
%     close(vidObj);

%% Localized admittance plot

    N0=(L+1)*(circle_rad/180);
    Sgt_local_mean=mean(Sgt_local,2);
    Stt_local_mean=mean(Stt_local,2);
    Sgg_local_mean=mean(Sgg_local,2);
    
    Z_local_mean=Sgt_local_mean./Stt_local_mean;
    Z_local_std=std(Sgt_local./Stt_local,0,2);
    Z_local_mean_std=Z_local_std/sqrt(NTapers);
    Degree=(1:numel(Z_local_mean))';
    Z_t(:,j)=Z_local_mean(GoodDegrees);
    Z_t_std(:,j)=Z_local_mean_std(GoodDegrees);
    
    Z_tm(:,(j-1)*NTapers+1:j*NTapers)=Sgt_local(GoodDegrees,:)./Stt_local(GoodDegrees,:);
    C_tm(:,(j-1)*NTapers+1:j*NTapers)=Sgt_local(GoodDegrees,:)./...
        sqrt(Stt_local(GoodDegrees,:).*Sgg_local(GoodDegrees,:));   
    
    C_local_mean=Sgt_local_mean./sqrt(Stt_local_mean.*Sgg_local_mean);
    C_local_std=std(Sgt_local./sqrt(Stt_local.*Sgg_local),0,2);
    C_local_mean_std=C_local_std/sqrt(NTapers);
    C_t(:,j)=C_local_mean(GoodDegrees);
    C_t_std(:,j)=C_local_mean_std(GoodDegrees);
 
    progressbar(j/numel(fii));
    j/numel(fii)

end
progressbar(1);