% close all
ccc;
tic;

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];


PhysicalConstants;

%% Input parameters
GM         = 62.8981e9;
Rref       = 476000;

a_obs      = 482.1130;
c_obs      = 446.1255;
J2_obs     = 0.0114;
T          = 9.074;

%% MCMC run parameters

N          = 10000; % number of iterations in a chain
Nc         = 12; % number of chains
s          = [1000, 50]; % step parameters

min_core_rad = 10000;
max_core_rad = 450000;
min_core_rho = 2160;
max_core_rho = 5000;

%% Basic calculations
M=GM/G;
f_obs=(a_obs-c_obs)/a_obs;
R=(a_obs.*a_obs.*c_obs).^(1/3);
V=1e9*4/3*pi*R.^3;

% sigma_fouter_obs=sqrt(...
%     (c_obs.*c_obs.*sigmaa.*sigmaa+a_obs.*a_obs.*sigmac.*sigmac)./(c_obs.^4));

router=(V./(4/3*pi)).^(1/3);

%% Make grid of values

ngrid   = 50;
rcore   = linspace(min_core_rad,max_core_rad,ngrid);
rhocoreg= linspace(min_core_rho,max_core_rho,ngrid);

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);
rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./...
    (4*pi*(rcorei.^3)-4*pi*(router^3));

%% Make initial values

data_obs=[f_obs J2_obs];

rhoouter_init = zeros(Nc,1)-1;
rcore_init    = zeros(Nc,1);
rhocore_init  = zeros(Nc,1);

for j=1:Nc
    while rhoouter_init(j) < 0        
        rcore_init(j)=rand*(max_core_rad-min_core_rad)+min_core_rad;
        rhocore_init(j)=rand*(max_core_rho-min_core_rho)+min_core_rho;        
        rhoouter_init(j)=-(3*M-4*pi*(rcore_init(j).^3).*rhocore_init(j))./...
            (4*pi*(rcore_init(j).^3)-4*pi*(router^3));        
    end 
end

param0=[rcore_init rhocore_init];
param=cell(1,Nc);

%% Run MCMC
tic
parfor j=1:Nc
    param_start = param0(j,:);
    param{j}    = mcmc_2l_run(N, s, data_obs, param_start, ...
        @cost_fun_flat,@step_param_flat,T,router,M);    
end
toc

param_all=[];

for i=1:Nc    
    param_all=[param_all; param{i}];    
end

%% Plot heatmap
% numbins = 100;
% marker = '.';
% markersize = 20;
% 
% outfile = heatscatter(param_all(:,1)/1000,param_all(:,2), ...
%     [], 'CeresHeatScatter.png', numbins, ...
%     markersize, marker, 1, 0, ...
%     'Core radius [km]', 'Core density [kg/m^{3}]', []);
% 
% xlim([min_core_rad max_core_rad]/1000);
% ylim([min_core_rho max_core_rho]);
% 
% set(gca,'FontSize',12);
% box on
% hold on



%% Plot another heat map

box_cond = (param_all(:,1) < min_core_rad) | ...
           (param_all(:,1) > max_core_rad) | ...
           (param_all(:,2) < min_core_rho) | ...
           (param_all(:,2) > max_core_rho);
       
param_all(box_cond,:) = [];       

n=100;
xi = linspace(min_core_rad,max_core_rad,n);
yi = linspace(min_core_rho,max_core_rho,n);

% xi = linspace(min(param_all(:,1)),max(param_all(:,1)),n);
% yi = linspace(min(param_all(:,2)),max(param_all(:,2)),n);

xr = interp1(xi,1:numel(xi),param_all(:,1),'nearest')';
yr = interp1(yi,1:numel(yi),param_all(:,2),'nearest')';

z = accumarray([xr' yr'], 1, [n n]);

[xii,yii]=meshgrid(xi,yi);

figure
xlim([min_core_rad max_core_rad]/1000);
ylim([min_core_rho max_core_rho]);
set(gca,'FontSize',12);
box on
hold on

pcolor(xii/1000,yii,z'); shading interp;
% plot(param_all(:,1)/1000,param_all(:,2),'.k');

levels = [1000 1250 1500 1750 2000 2500 3000];
[C,h] = contour(rcorei/1000,rhocorei,rhoouteri,levels,...
    'Color',[1 1 0],'LineWidth',2); shading interp;

text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.4 .4 .4]);

%%

rhoouter_mc=-(3*M-4*pi*(param_all(:,1).^3).*param_all(:,2))./...
    (4*pi*(param_all(:,1).^3)-4*pi*(router^3));

n=200;

st = (router-param_all(:,1))/1000;

xi = linspace(min(rhoouter_mc),max(rhoouter_mc),n);
yi = linspace(min(st),max(st),n);

xr = interp1(xi,1:numel(xi),rhoouter_mc,'nearest')';
yr = interp1(yi,1:numel(yi),st,'nearest')';

z = accumarray([xr' yr'], 1, [n n]);

[xii,yii]=meshgrid(xi,yi);

fig_shell=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

pcolor(xii,yii,(z')); shading interp;
% plot(rhoouter_mc,st,'.k');

xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);
xlim([min(rhoouter_mc) max(rhoouter_mc)]);
ylim([min(st),max(st)]);











