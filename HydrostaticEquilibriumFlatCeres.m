
% ccc;
tic;
%% Input parameters


GM=62.8981e9;

G=6.67e-11;
M=GM/G;

% V=4.5228e8*1e9;

Rref=500000;

% Equatorial semi-axis (km)   487.30  1.8
% Polar semi-axis (km)        454.70  1.6


a=487.3;
c=454.7;
sigmaa=1.8;
sigmac=1.6;

% a=479.7;
% c=444.4;
% sigmaa=2.3;
% sigmac=2.1;

R=(a*a*c)^(1/3);

V=1e9*4/3*pi*R^3;

fouter_obs=(a-c)/a;
sigma_fouter_obs=sqrt(...
    (c*c*sigmaa*sigmaa+a*a*sigmac*sigmac)./(c.^4));


T=9.074;

router=(V/(4/3*pi))^(1/3);

% T=9:0.1:10;


% W=2*pi/(T*3600);

rcore=2500:5000:450000;

fouter_obs=0.066899240714139;
sigma_fouter_obs=0.0047;

fcore0=0.1;
fouter0=0.1;

% rhocore=4000:25:8000;

%% Plotting settings
fig1=figure; hold on; 
ax1=gca;
set(ax1,'FontSize',12);

data=load('ListJ2.txt');

xlim([10 450]);
ylim([2000 6000]);

xlabel('Core size [km]','FontSize',12);
ylabel('Core density [kg/m^3]','FontSize',12);

rhocore0=2000;
box on;

%% Grid of core radii and densities
rhocoreg=2000:100:6000;

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);

rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));


progressbar(0)

for i=1:numel(rhocorei)

    [fhi(i,:),fval]=HydrostaticStateExact2l(router,rcorei(i),T,rhoouteri(i),rhocorei(i),0.1,0.1);
    fouter0=fhi(i,1);
    fcore0=fhi(i,2);
    
    progressbar(i/numel(rhocorei));
end

progressbar(1);

fcorei=fhi(:,2);
fouteri=fhi(:,1);

fcorei=reshape(fcorei,size(rhocorei));
fouteri=reshape(fouteri,size(rhocorei));

% surf(rcorei/1000,rhocorei,fcorei); shading interp
% caxis([0.08 0.16]);

% alpha(0.1);
%% Contour mantle density
levels=[0:500:1500 1750:250:4000];

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri,levels,'Color',[0.4 0.4 0.4]);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7]);

box on;

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri,[1000 1000],'Color',[0.4 0.4 0.4],'LineWidth',3);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7]);

pcolor(rcorei/1000,rhocorei,rhoouteri); shading interp;
cbar=colorbar('FontSize',12);
ylabel(cbar,'Mantle density [kg/m^3] ','FontSize',12);
caxis([0 2000]);

%% Contour core flattening

% levels=0:0.01:0.6;
% [C,h]=contour(rcorei/1000,rhocorei,fcorei,levels,'Color',[0 1 0]);
% text_handle = clabel(C,h);


%% Contour outerflattening

levels=0:0.005:0.6;
[C,h]=contour(rcorei/1000,rhocorei,fouteri,levels,'Color',[1 1 1],'LineWidth',2);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 1],...
    'Edgecolor',[.7 .7 .7]);

levels=[fouter_obs-sigma_fouter_obs fouter_obs+sigma_fouter_obs];
 
[C,h]=contour(rcorei/1000,rhocorei,fouteri,levels,'Color','k','LineWidth',4);
% text_handle = clabel(C,h);

% p = findobj(h, 'cdata', max(levels));
% set(p, 'facecolor', 'white');

%% Computing hydrostatic J2

J2hi=RadFlat2J2(router,rcorei,fouteri,fcorei,rhoouteri,rhocorei,Rref);

nmax=10;
% 
% J=HydroCoeffs2Layer(router,rcorei,fouteri,fcorei,rhoouteri,rhocorei,Rref,nmax);

figure; hold on;

levels=0:2.5e-4:0.004;
% [C,h]=contour(rcorei/1000,rhocorei,J2hi,levels,'Color',[0 0 1]);
[C,h]=contour(rcorei/1000,rhocorei,J2hi,'Color',[0 1 1],'LineStyle','-','LineWidth',2);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 1],...
    'Edgecolor',[.7 .7 .7]);

% J2exp_min=griddata(rcorei/1000,rhocorei,J(:,:,1),C(1,:),C(2,:),'cubic');
% J2exp_max=griddata(rcorei/1000,rhocorei,J(:,:,1),C(1,:),C(2,:),'cubic');

% [C,h]=contour(rcorei/1000,rhocorei,J2hi,'Color',[0 1 1],'LineStyle','-','LineWidth',2);


% nc=50;
% mc=30;
% 
% for i=1:size(J,3)
%     Jslice(i)=J(nc,mc,i);
% end
% 
% figure; hold on;
% set(gca,'FontSize',20,'YScale','log');
% 
% 
% plot(2:2:nmax,Jslice,'-o','LineWidth',3,'MarkerSize',6);
% 
% xlabel('Degree','FontSize',20);
% ylabel('J_{2n} ','FontSize',20);
% 
% 
% for i=1:size(rhocorei,1)
%     for j=1:size(rhocorei,2)
%         
%         for k=1:size(J,3)
%             Jslice(k)=J(i,j,k);
%         end
%         
%         p=polyfit(2:2:nmax,Jslice,1);
%         
%         pe(i,j)=p(1);
%         
%     end   
% end


% 
% figure; hold on;
% set(gca,'FontSize',20);
% pcolor(rcorei/1000,rhocorei,J(:,:,1)); shading interp;
% colorbar('FontSize',20);
% 
% 
% xlabel('Core size [km] ','FontSize',12);
% ylabel('Core density [kg/m^3] ','FontSize',12);




%% Computing C/MR^2
% 
% [a1,c1]=f2axes(router,fouteri);
% [a2,c2]=f2axes(rcorei,fcorei);
% 
% M1=4/3*pi*rhoouteri*router.^3;
% M2=4/3*pi*(rhocorei-rhoouteri).*rcorei.^3;
% 
% M=M1+M2;
% 
% Ch1=0.2*(M1).*(a1.^2+c1.^2);
% Ch2=0.2*(M2).*(a2.^2+c2.^2);
% 
% Ch=Ch1+Ch2;
% 
% lambdah=Ch./(M.*router.^2);
% 
% levels=-0.5:0.005:0.6666;
% [C,h]=contour(rcorei/1000,rhocorei,lambdah,levels,'Color',[0 0 0],'LineWidth',2);
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 1],...
%     'Edgecolor',[.7 .7 .7]);


% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7]);

%% Plotting J2 solution

% rcore_m=85:1:230;
% 
% plot(rcore_m,data(:,1),'r')
% plot(rcore_m,data(:,2),'r')


% f1obs=0.19;
% f2max=0.163;
% f2min=0.098;
% 
% rhocore_j2soln_1=CoreDensityFromJ2Mass(M,rcore,f2max,Rref,J2outer,J2obs,V);
% rhocore_j2soln_2=CoreDensityFromJ2Mass(M,rcore,f2min,Rref,J2outer,J2obs,V);

% 


% plot(ax1,rcore/1000,rhocore_j2soln_1,'r')
% plot(ax1,rcore/1000,rhocore_j2soln_2,'r')
% 
% X=[rcore/1000,fliplr(rcore/1000)];                %#create continuous x value array for plotting
% Y=[rhocore_j2soln_1 fliplr(rhocore_j2soln_2)];              %#create y values for out and then back
% fill(X,Y,'r'); 

% %% Main loop
% 
% Tol=0.0001;
% Eps=Tol+1;
% 
% progressbar(0, 0)
% 
% for k=1:numel(T)
%     
% 
% fcore0=0.1;
% fouter0=0.1;
% 
% 
%     for i=1:numel(rcore)
%         
%    
%     
%     rhocore_h(i)=fzero(@(rhocore) FunFlat(rhocore,M,rcore(i),router,T(k),fouter_obs),rhocore0);
%     
%     rhoouter=-(3*M-4*pi*((rcore(i)).^3).*rhocore_h(i))./(4*pi*(rcore(i).^3)-4*pi*(router^3));
%     
%     [fh,d2]=HydrostaticStateExact2l(router,rcore(i),T(k),rhoouter,rhocore_h(i),fouter0,fcore0);
%     
%     rhocore_h2(i)=CoreDensityFromJ2Mass(M,rcore(i),fh(2),Rref,J2outer,J2obs,V);
%     
%     
%         while (Eps>Tol)
%             
%         rhoouter=-(3*M-4*pi*((rcore(i)).^3).*rhocore_h2(i))./(4*pi*(rcore(i).^3)-4*pi*(router^3));
%             
%         fhold=fh;
%         
%         [fh,d2]=HydrostaticStateExact2l(router,rcore(i),T(k),rhoouter,rhocore_h2(i),fh(1),fh(2));
%         
%         rhocore_h2(i)=CoreDensityFromJ2Mass(M,rcore(i),fh(2),Rref,J2outer,J2obs,V);
%         
%         Eps=max(abs(fhold-fh));
%         
%         end
%      
%     rhocore0=rhocore_h(i);
%     progressbar(k/numel(T), i/numel(rcore));
%     
%     fouter0=fh(1);
%     fcore0=fh(2);
%     
%     
%     end
%     
%     
% 
%     
%     % plot solution for given T
%     
%     plot(ax1,rcore/1000,rhocore_h,'-k','LineWidth',1);
%     plot(ax1,rcore/1000,rhocore_h2,'-b','LineWidth',1);
%     
%     progressbar(k/numel(T), i/numel(rcore));
%     
%     rhocore0=rhocore_h(1);
% 
% 
% end
% 
% progressbar(1,1);

%% Plotting




% set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
% set(gcf, 'PaperPositionMode','auto')
% 
% box on
% 
% 
% print(fig1, '-dpsc', 'HydrostaticEquilibriumExactSolutionJ2flat.eps');


toc







