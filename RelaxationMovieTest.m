ccc

%% Image settings
fntsize = 20;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

ccj = {[0.0 0.3 0.9],[0.8 0.1 0.3]};

%% Body parameters
r   = 470;
rho = 2161;
T   = 7;
L   = 40;
Rref = 470;

%% Read data

movie_filename ='RelaxationMovie_large_contrast_4.avi';
folder_path   = '/Users/antonermakov/Dawn/FE/output/output_1';
filename_mesh = getAllFiles(folder_path,'_mesh');
filename_surf = getAllFiles(folder_path,'_surface');

data = load([folder_path '/physical_times.txt']);
t = data(:,2)/(365.2422*86400);

%% Hydrostatic equilibrium computation
[fh,fval]=HydrostaticStateExact(r*1000,T,rho,0.1);
[fh,fval]=HydrostaticStateExact2l(r*1000,410000,T,1465,2491,0.1, 0.1);

[a,c]=f2axes(r,fh(1));

ang = linspace(0,pi/2,100);
xell = a*cos(ang);
zell = c*sin(ang);

%% Figure for the movie


relax_pl = figure('Position',[1 1 1200 500]); 
pl_shape = subplot(1,2,1);
hold on;

xlim([0 500]);
ylim([0 500]);
xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
ylabel('z [km]','FontSize',fntsize,'interpreter','latex');
box on;

pl_spectrum = subplot(1,2,2);
hold on; box on; grid on;
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'FontSize',fntsize);
xlabel('Frequency [cycles/km]','FontSize',fntsize,'interpreter','latex');
ylabel('Topography power [$\textrm{km}^{3}$]','FontSize',fntsize,...
    'interpreter','latex');

ylim([1e2 1e8]);
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
plot(xell,zell,'r--','LineWidth',4);

xlim([0 500]);
ylim([0 500]);

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

plot(xell,zell,'r--','LineWidth',4);


subplot(pl_spectrum);

lmcosi_limb = quad2plm(filename_surf{2},L);
[sdl_limb,l_limb] = plm2spec(lmcosi_limb);

lambda=2*pi./l_limb;
lambda_linear=lambda*Rref;
kf=1./lambda_linear;

plot_ceres = plot(kf(1:2:end),sdl_limb(1:2:end),...
    '-o','MarkerSize',2,'Color','k');

plot(kf(1:2:end),sdl_limb(1:2:end),...
    '-o','MarkerSize',2,'Color','r');

v = VideoWriter(movie_filename);
open(v);

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
    
    lmcosi_limb = quad2plm(filename_surf{i},L);
    [sdl_limb,l_limb] = plm2spec(lmcosi_limb);
    set(plot_ceres,'YData',sdl_limb(1:2:end));
       
    frame = getframe(gcf);
    writeVideo(v,frame);
    i/numel(filename_surf)
    
end

close all
close(v)



