ccc

%% Image settings
fntsize = 20;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

ccj = {[1 0 0 ],...
       [0 1 0],...
       [0 0 1],...
       [0.6 0.6 0.6],...
       [0.6 0.6 0.6],...
       [0.6 0.6 0.6]};

%% Body parameters
r    = 470;
rho  = 2200;
T    = 90000.0;
L    = 80;
Rref = 470;

%% Read data

movie_filename ='RelaxationMovie_V2_new_ell_inv_eta.avi';
folder_path   = '/Users/antonermakov/Dawn/FE/output/output_1';
filename_mesh = getAllFiles(folder_path,'_mesh');
filename_surf = getAllFiles(folder_path,'_surface');
filename_pl   = getAllFiles(folder_path,'_failurelocations00');
filename_ps   = getAllFiles(folder_path,'_principalstresses00');
filename_viscbase  = getAllFiles(folder_path,'_baseviscosities');
filename_viscreg   = getAllFiles(folder_path,'_viscositiesreg00');
filename_stress    = getAllFiles(folder_path,'_stresstensor00');
filename_flow      = getAllFiles(folder_path,'_flow00');

data = load([folder_path '/physical_times.txt']);
t = data(:,2);

%% Hydrostatic equilibrium computation
% [fh,fval]=HydrostaticStateExact(r*1000,T,rho,0.1);
[fh,fval]=HydrostaticStateExact2l(r*1000,r*1000-150000,T,2200,2200,0.1, 0.1);

[a,c]=f2axes(r,fh(1));
[a_cmb,c_cmb]=f2axes(r-150,fh(2));

ang = linspace(0,pi/2,100);
xell = a*cos(ang);
zell = c*sin(ang);

xell_cmb = a_cmb*cos(ang);
zell_cmb = c_cmb*sin(ang);

%% Figure for the movie

relax_pl = figure('Position',[1 1 1.2*1200 1.2*500],'Color','w');
pl_shape = subplot(1,2,1);
hold on;

xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
ylabel('z [km]','FontSize',fntsize,'interpreter','latex');
box on;

pl_spectrum = subplot(1,2,2);
hold on; box on; grid on;
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'FontSize',fntsize);
xlabel('Frequency [cycles/km]','FontSize',fntsize,'interpreter','latex');
ylabel('Spectral density [$\textrm{km}^{2}$]','FontSize',fntsize,...
    'interpreter','latex');

% load and plot read Ceres spectral density
load('Spectral_Density_Ceres.mat');

lambda_Ceres=2*pi./l_Ceres;
lambda_linear_Ceres=lambda_Ceres*Rref;
kf_Ceres=1./lambda_linear_Ceres;

plot(kf_Ceres,sdl_Ceres,'-b','LineWidth',3);

ylim([1e0 1e8]);
xlim([1e-4 1e-1]);

subplot(pl_shape)
% fig_spec=figure('Color','w','Position',[1 1 1000 1000]);
% set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on; axis equal;

xlabel('r [km]','FontSize',fntsize,'interpreter','latex');
ylabel('z [km]','FontSize',fntsize,...
    'interpreter','latex');

set(gca, 'FontSize',fntsize);
set(gca,'XTick',0:100:500);
set(gca,'YTick',0:100:500);

title(['t = ' num2str(t(2),'%6.2e') ' [y]'],'FontSize',fntsize,'interpreter','latex');

xlim([0 600]);
ylim([0 600]);

meshStruct = Read_ucd(filename_mesh{2});
V = meshStruct.V/1000;
E = meshStruct.E;
cell_mat = meshStruct.cell_mat;

% draw and record first frame
for j=1:size(E,1)
    
    color = ccj{cell_mat(j)+1};
    p(j) = patch([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],color);
    l(j) = line([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],'Color','k');
end

plot(xell,zell,'k--','LineWidth',4);
plot(xell_cmb,zell_cmb,'k--','LineWidth',4);

% base viscocities
bv = load(filename_viscbase{2});
x = bv(:,1);
z = bv(:,2);

% plot faillure location
fl = load(filename_pl{2});
try
    failure_plot = plot(fl(:,1)/1000,fl(:,2)/1000,'xy','MarkerSize',10);
end
% plot stresses
s = load(filename_ps{2});
% p_stresses_plot = scatter(x/1000,z/1000,10,s(:,1)./s(:,2),'filled');
% caxis([0.2 5]);

subplot(pl_spectrum);

lmcosi_limb = quad2plm(filename_surf{2},L);
[sdl_limb,l_limb] = plm2spec(lmcosi_limb);

lambda=2*pi./l_limb;
lambda_linear=lambda*Rref;
kf=1./lambda_linear;

plot_ceres = plot(kf(1:2:end),sdl_limb(1:2:end),...
    '-o','MarkerSize',2,'Color','k','LineWidth',3);

plot(kf(1:2:end),sdl_limb(1:2:end),...
    '-o','MarkerSize',2,'Color','r','LineWidth',3);

v = VideoWriter(movie_filename);
open(v);

legend({'Real Ceres','Relaxed','Power Law'},'FontSize',fntsize);

frame = getframe(gcf);
writeVideo(v,frame);

% draw and record all other frames
for i=3:1:numel(filename_mesh)-1
    
    meshStruct = Read_ucd(filename_mesh{i});
    V = meshStruct.V/1000;
    E = meshStruct.E;
    
    subplot(pl_shape);
    title(['t = ' num2str(t(i),'%6.2e') ' [y]'],'FontSize',fntsize,'interpreter','latex');
    
    for j=1:size(E,1)
        set(l(j),'XData',[V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)]);
        set(l(j),'YData',[V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
        set(p(j),'XData',[V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)]);
        set(p(j),'YData',[V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
    end
    
    try
        fl = load(filename_pl{i});
        set(failure_plot,'XData',fl(:,1)/1000,'YData',fl(:,2)/1000);
    end
    
    lmcosi_limb = quad2plm(filename_surf{i},L);
    [sdl_limb,l_limb] = plm2spec(lmcosi_limb);
    set(plot_ceres,'YData',sdl_limb(1:2:end));
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    i/numel(filename_surf)
    
end

close all
close(v)



