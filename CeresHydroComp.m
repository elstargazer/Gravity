% close all
ccc;
tic;
%% Input parameters


GM=62.8981e9;

G=6.67e-11;
M=GM/G;

% V=4.5228e8*1e9;

Rref=500000;

% Equatorial semi-axis (km)   487.30  1.8
% Polar semi-axis (km)        454.70  1.6


a=[487.3 479.7 483.4990 482.1130];
c=[454.7 444.4 447.6 446.1255];
sigmaa=[1.8 2.3 0.7 0];
sigmac=[1.6 2.1 0.4 0];

% a=;
% c=;
% sigmaa=2.3;
% sigmac=2.1;

R=(a.*a.*c).^(1/3);

V=1e9*4/3*pi*R.^3;

fouter_obs=(a-c)./a;
sigma_fouter_obs=sqrt(...
    (c.*c.*sigmaa.*sigmaa+a.*a.*sigmac.*sigmac)./(c.^4));


T=9.074;

router=(V./(4/3*pi)).^(1/3);

router=router(4);

rhomean=M/V(4)
% T=9:0.1:10;


% W=2*pi/(T*3600);

fcore0=0.1;
fouter0=0.1;

% f_d=[0.0598 0.0651]

% rhocore=4000:25:8000;


%% Grid of core radii and densities
rcore=2500:2500:470000;
rhocoreg=2000:100:6000;

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);

rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));

rhoouteri(rhoouteri<0)=NaN;


progressbar(0)

for i=1:numel(rhocorei)

    [fhi(i,:),fval]=HydrostaticStateExact2l(router,rcorei(i),T,rhoouteri(i),rhocorei(i),0.1,0.1);
%     fhi_an(:,i)=HydrostaticState2LayerAn(rcorei(i),router,T,rhocorei(i),rhoouteri(i),0.1,0.1,'RD');
    fouter0=fhi(i,1);
    fcore0=fhi(i,2);
    
    progressbar(i/numel(rhocorei));
end

progressbar(1);

fcorei=fhi(:,2);
fouteri=fhi(:,1);

fcorei=reshape(fcorei,size(rhocorei));
fouteri=reshape(fouteri,size(rhocorei));

% fouteri_an=reshape(fhi_an(2,:),size(rhocorei));

% surf(rcorei/1000,rhocorei,fcorei); shading interp
% caxis([0.08 0.16]);

% alpha(0.1);
%% Plotting settings
fig1=figure; hold on; 
ax1=gca;
set(ax1,'FontSize',12);

%data=load('ListJ2.txt');

xlim([10 450]);
ylim([2000 6000]);

xlabel('Core size [km]','FontSize',20);
ylabel('Core density [kg/m^3]','FontSize',20);

rhocore0=2000;
box on;

%% Contour mantle density
levels=[0:500:1500 1750:250:4000];

% [C,h]=contour(rcorei/1000,rhocorei,rhoouteri,levels,'Color',[0.4 0.4 0.4]);
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7]);

box on;

% [C,h]=contour(rcorei/1000,rhocorei,rhoouteri,[1000 1000],...
%     'Color',[0.4 0.4 0.4],'LineWidth',3);
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7]);

pcolor(rcorei/1000,rhocorei,rhoouteri); shading interp;
cbar=colorbar('FontSize',20);
ylabel(cbar,'Mantle density [kg/m^3] ','FontSize',20);
caxis([0 2000]);

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri,[1000 1250 1500 1750 2000 2500 3000],...
    'Color',[1 1 0],'LineWidth',2); shading interp;
cbar=colorbar('FontSize',12);
ylabel(cbar,'Mantle density [kg/m^3] ','FontSize',12);

text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.4 .4 .4]);


caxis([0 2000]);

%% Contour core flattening

% levels=0:0.01:0.6;
% [C,h]=contour(rcorei/1000,rhocorei,fcorei,levels,'Color',[0 1 0]);
% text_handle = clabel(C,h);


%% Contour outerflattening

levels=0:0.005:0.6;
[C,h]=contour(rcorei/1000,rhocorei,fouteri,levels,'Color',[0.85 0.85 0.85],'LineWidth',2);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[0.85 0.85 0.85],...
    'Edgecolor',[.7 .7 .7]);



rho_s=3300;
rho_i=920;




rcorei
rhocorei
rhoouteri


x1=(rho_s-rhocorei)./(rho_s-rho_i); % ice in core
x2=(rho_s-rhoouteri)./(rho_s-rho_i); % ice in mantle

rd3=(router.^3-rcorei.^3);

R=(x1.*rho_i.*rcorei.^3+x2.*rho_i.*rd3)./...
    ((1-x1).*rho_s.*rcorei.^3+(1-x2).*rho_s.*rd3);


[C,h]=contour(rcorei/1000,rhocorei,R,'Color',[0.55 0.85 0.45],'LineWidth',2);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[0.85 0.85 0.85],...
    'Edgecolor',[.7 .7 .7]);
    


fig1=figure; hold on; 
ax1=gca;
set(ax1,'FontSize',12);


pcolor(rcorei/1000,rhocorei,R); shading interp

caxis([0 1]);

colorbar;

xlim([10 450]);
ylim([2000 6000]);



(rho_s-rhomean)./(rho_s-rho_i)






















