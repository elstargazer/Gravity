% ccc

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

folder_path = '/Users/antonermakov/Dawn/FE/output/output_1';
filename_quad = getAllFiles(folder_path,'_mesh');

fig_spec=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
ylabel('y [km]','FontSize',fntsize,...
    'interpreter','latex');

v = VideoWriter('RelaxationMovie.avi');
open(v);

data = load([folder_path '/physical_times.txt']);
t = data(:,2)/(365.2422*86400);

for i=1:numel(filename_quad);
    
    meshStruct = Read_ucd(filename_quad{i});
    V = meshStruct.V/1000;
    E = meshStruct.E;
    
    xlim([0 500]);
    ylim([0 500]);
    xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
    ylabel('y [km]','FontSize',fntsize,...
        'interpreter','latex');
    grid on; box on;
    
    title(['t = ' num2str(t(i),'%6.2e')],'FontSize',fntsize,'interpreter','latex');
    
    
    for j=1:size(E,1)
        
        line([V(E(j,1),1) V(E(j,2),1)], [V(E(j,1),2) V(E(j,2),2)],'Color','k');
        line([V(E(j,2),1) V(E(j,3),1)], [V(E(j,2),2) V(E(j,3),2)],'Color','k');
        line([V(E(j,3),1) V(E(j,4),1)], [V(E(j,3),2) V(E(j,4),2)],'Color','k');
        line([V(E(j,4),1) V(E(j,1),1)], [V(E(j,4),2) V(E(j,1),2)],'Color','k');
        
    end
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    clf;
    
    i/numel(filename_quad)
    
end

close(v)
