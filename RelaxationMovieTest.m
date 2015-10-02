ccc

%% Image settings
fntsize = 20;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

%% Body parameters
r   = 470;
rho = 2161;
T   = 7;

%% Read data
folder_path = '/Users/antonermakov/Dawn/FE/output/output_1';
filename_quad = getAllFiles(folder_path,'_mesh');

data = load([folder_path '/physical_times.txt']);
t = data(:,2)/(365.2422*86400);

%% Hydrostatic equilibrium computation
[fh,fval]=HydrostaticStateExact(r*1000,T,rho,0.1);
[a,c]=f2axes(r,fh);

ang = linspace(0,pi/2,100);
xell = a*cos(ang);
zell = c*sin(ang);

%% Figure for the movie
fig_spec=figure('Color','w','Position',[1 1 1000 1000]);
% set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on; axis equal;

xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
ylabel('y [km]','FontSize',fntsize,...
    'interpreter','latex');

v = VideoWriter('RelaxationMovie_6.avi');
open(v);


for i=2:numel(filename_quad)-1
    
    meshStruct = Read_ucd(filename_quad{i});
    V = meshStruct.V/1000;
    E = meshStruct.E;
    

    xlabel('r [km]','FontSize',fntsize,'interpreter','latex');
    ylabel('z [km]','FontSize',fntsize,...
        'interpreter','latex');
    set(gca, 'FontSize',fntsize);
    set(gca,'XTick',0:100:500);
    set(gca,'YTick',0:100:500);

    grid on; box on;
    axis equal;
    
    xlim([0 500]);
    ylim([0 500]);
    
    title(['t = ' num2str(t(i),'%6.2e') ' [y]'],'FontSize',fntsize,'interpreter','latex');
        plot(xell,zell,'r--','LineWidth',4);

    
    for j=1:size(E,1)
        
        line([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
             [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],'Color','k');
%         line([V(E(j,2),1) V(E(j,3),1)], [V(E(j,2),2) V(E(j,3),2)],'Color','k');
%         line([V(E(j,3),1) V(E(j,4),1)], [V(E(j,3),2) V(E(j,4),2)],'Color','k');
%         line([V(E(j,4),1) V(E(j,1),1)], [V(E(j,4),2) V(E(j,1),2)],'Color','k');
        
    end
    

    frame = getframe(gcf);
    writeVideo(v,frame);
    
    cla;
    
    i/numel(filename_quad)
    
end

close(v)



