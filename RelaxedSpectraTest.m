close all;

%% Plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

%% reading data
% filename_quad = {  '../CeresFE/FE/outer_vertices_0.txt',...
%     '../CeresFE/FE/outer_vertices_1.txt',...
%      '../CeresFE/FE/outer_vertices_2.txt'};

folder_path = '../CeresFE/FE/Output/';
 
filename_quad = getAllFiles(folder_path);

L = 45;

fig_spec=figure; 
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('Degree','FontSize',fntsize);
ylabel('Power [km^2]','FontSize',fntsize);

ylim([1e2 1e8]);
xlim([1 L]);


ccj = cool(numel(filename_quad));

% figure; hold on;
% xlim([0 pi/2]);
% 
% for i = 1:numel(filename_quad)  
%     data = load(filename_quad{i}); 
%     
%     [~,lat,r] = cart2sph(data(:,1),0,data(:,2));
%     
% %     plot(data(:,1),data(:,2),'.','Color',ccj(i,:));
% %     plot(data(:,2),data(:,1),'.','Color',ccj(i,:));
% 
%     plot(lat,r,'.','Color',ccj(i,:));
% end
% 
% ylim([450000 485000]);
for i = 1:numel(filename_quad)
    lmcosi_limb = quad2plm(filename_quad{i},L);
    [sdl_limb,l_limb] = plm2spec(lmcosi_limb);
    figure(fig_spec);
    pl_relax = plot(l_limb(1:2:end),sdl_limb(1:2:end)...
        ,'-o','MarkerSize',2,'Color',ccj(i,:));
    
    drawnow;
end





