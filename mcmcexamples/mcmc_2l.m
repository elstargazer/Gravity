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

a_obs=482.1130;
c_obs=446.1255;

f_obs=(a_obs-c_obs)/a_obs;
J2_obs=0.01163;

% a=;
% c=;
% sigmaa=2.3;
% sigmac=2.1;

R=(a_obs.*a_obs.*c_obs).^(1/3);

V=1e9*4/3*pi*R.^3;


% sigma_fouter_obs=sqrt(...
%     (c_obs.*c_obs.*sigmaa.*sigmaa+a_obs.*a_obs.*sigmac.*sigmac)./(c_obs.^4));

T=9.074;

router=(V./(4/3*pi)).^(1/3);

rcore=2500:2500:470000;
rhocoreg=2000:100:6000;

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);
rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));
rhoouteri(rhoouteri<0)=NaN;

%% 

N=50000; % number of iterations in a chain
Nc=12; % number of chains
s=[1000, 50]; % step parameters


min_core_rad=10000;
max_core_rad=450000;
min_core_rho=2200;
max_core_rho=5000;

data_obs=[f_obs J2_obs];

HowMuch=10;
rcore_init=rand(HowMuch*Nc,1)*(max_core_rad-min_core_rad)+min_core_rad;
rhocore_init=rand(HowMuch*Nc,1)*(max_core_rho-min_core_rho)+min_core_rho;

rhoouter_init=-(3*M-4*pi*(rcore_init.^3).*rhocore_init)./...
    (4*pi*(rcore_init.^3)-4*pi*(router^3));

Norig=numel(rcore_init);
for i=1:Norig
    try
        if (rhoouter_init(i)<0)
            rcore_init(i)=[];
            rhocore_init(i)=[];
            rhoouter_init(i)=[];
        end
    catch
        
    end
end

rcore_init=rcore_init(1:Nc);
rhocore_init=rhocore_init(1:Nc);
rhoouter_init=rhoouter_init(1:Nc)

param0=[rcore_init rhocore_init];
param=cell(1,Nc);

tic
parfor j=1:Nc
    param_start = param0(j,:);
    param{j}=mcmc_2l_run(N, s, data_obs, param_start, ...
        @cost_fun_flat,@step_param_flat,T ,router ,M);    
end
toc

param_all=[];

for i=1:Nc    
    param_all=[param_all; param{i}];    
end

% figure; hold on;
% plot(param_all(:,1)/1000,param_all(:,2),'.')
% xlim([0 470]);
% ylim([2000 6000]);


%% Plot Heatmap
numbins = 100;
marker = 'o';
markersize = 10;
outfile = heatscatter(param_all(:,1)/1000,param_all(:,2), ...
    [], 'CeresHeatScatter.png', numbins, ...
    markersize, marker, 1, 0, ...
    'Core radius [km]', 'Core density [kg/m^{3}]', []);

xlim([0 470]);
ylim([2000 6000]);

set(gca,'FontSize',12);
box on
hold on

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri,[1000 1250 1500 1750 2000 2500 3000],...
    'Color',[1 1 0],'LineWidth',2); shading interp;

text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.4 .4 .4]);









