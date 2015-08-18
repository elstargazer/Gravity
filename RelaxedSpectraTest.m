ccc;

%% Plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';

%% Load Ceres shape

load CeresShapeSPCSurvey.mat;
lmcosi_ceres(4,3) = 0;
[sdl_ceres,l_ceres] = plm2spec(lmcosi_ceres);

%% reading data
% filename_quad = {  '../CeresFE/FE/outer_vertices_0.txt',...
%     '../CeresFE/FE/outer_vertices_1.txt',...
%      '../CeresFE/FE/outer_vertices_2.txt'};

folder_path = '~/Dawn/FE/output/';

L = 60;

% t = [0:1]*1e3;
% t = [0:8 10:19]*1e3;

filename_quad = getAllFiles(folder_path);

Ncolors = 128;
ccj = cool(Ncolors);

for i=1:numel(filename_quad)
    [~,time_str,~]=fileparts(filename_quad{i});
    time_str(5:6)
    t(i) = str2double(time_str(5:6));
end

min_t = min(t);
max_t = max(t);

fig_spec=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('Degree','FontSize',fntsize,'interpreter','latex');
ylabel('Power spectral density [$\textrm{m}^{2}$]','FontSize',fntsize,...
    'interpreter','latex');

ylim([1e-2 1e8]);
xlim([1 360]);

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
% ylim([450000 485000]);

lmcosi_limb_0 = quad2plm(filename_quad{1},L);
[sdl_limb_0,l_limb_0] = plm2spec(lmcosi_limb_0);

figure(fig_spec);
plot_ceres = plot(l_ceres(3:end),sdl_ceres(3:end),...
    '-o','MarkerSize',2,'Color','k');

plot_fit = plot(l_limb_0(1:2:end),sdl_limb_0(1:2:end),...
    '-o','MarkerSize',2,'Color','r');

for i = 2:numel(filename_quad)
    
    [~,time_str,~]=fileparts(filename_quad{i});
    t(i) = str2double(time_str(5:6));
    lmcosi_limb = quad2plm(filename_quad{i},L);
    [sdl_limb(:,i-1),l_limb] = plm2spec(lmcosi_limb);
    
    d = (log10(sdl_limb_0) - log10(sdl_limb(:,i-1)));
    d(3:2:end)
    
    tau(:,i-1) = t(i)./(log(sdl_limb_0) - log(sdl_limb(:,i-1)));
    
    figure(fig_spec);
    
    color_ind = fix((t(i)-min_t)/(max_t-min_t)*(Ncolors-1))+1;
    pl_relax = plot(l_limb(1:2:end),sdl_limb(1:2:end,i-1)...
        ,'-o','MarkerSize',2,'Color',ccj(color_ind,:));
       
%     pl_relax_homo = plot(l_limb(1:2:end),sdl_limb(1:2:end,i-1)...
%         ,'-o','MarkerSize',2,'Color',[0 0.7 0] );
    
    drawnow;
end

% PrintWhite(fig_spec,[fig_folder 'Fig_layer_relaxed_spectrum.jpg']);

%% Relxation time
% plot_relax_times = figure;
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on;grid on; box on;
% set(gca,'YScale','log');
% 
% ccj = cool(size(tau,2));
% 
% for i = 1:size(tau,2)
%     
%     plot(2:2:L,tau(3:2:end,i),'o-k','LineWidth',3,'MarkerSize',5,...
%         'Color',ccj(i,:));
%     
% end
% 
% log10(tau(3:2:end,:))
% 
% xlabel('Degree','FontSize',fntsize,'interpreter','latex');
% ylabel('Characteristic relaxation time [years]','FontSize',fntsize,...
%     'interpreter','latex');
% 
% % figure;
% % set(gcf, 'Units','centimeters', 'Position',im_size)
% % set(gcf, 'PaperPositionMode','auto')
% % set(gca, 'FontSize',fntsize);
% % hold on;grid on; box on;
% 
% tol_order = 1;
% 
% for ind = 3:2:size(sdl_limb,1)
%     
%     l_fit = l_limb(ind);
%     
%     
%     cond_order = ((log10(sdl_limb_0(ind)) - log10(sdl_limb(ind,:))) < tol_order);
%     first_false_ind = find(~cond_order,1);
%     cond_order(first_false_ind:end) = 0;
%     
%     %     [p,S] = polyfit(t(cond_order),...
%     %         [log(sdl_limb_0(ind,cond_order)) log(sdl_limb(ind,cond_order))],1);
%     [p,S] = polyfit([0 t(find(cond_order)+1)],[log(sdl_limb_0(ind,:)) ...
%         log(sdl_limb(ind,cond_order))],1);
%     
%     tau_fit(ind) = -1/p(1);
%     
%     %     cf = fit(t',[log(sdl_limb_0(ind,:)) log(sdl_limb(ind,:))]','poly1');
%     %
%     %     cf_coeff = coeffvalues(cf);
%     %     cf_confint = confint(cf);
%     %     a = cf_coeff(1);
%     %     b = cf_coeff(2);
%     %     a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
%     %     b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
%     %
%     %     tau_uncert = sqrt((a_uncert^2)/(a.^4));
%     
%     plot(ind-1,tau_fit(ind),'o','MarkerSize',6,...
%         'MarkerFaceColor','r','MarkerEdgeColor','k');
%     %     errorbar(ind-1,tau_fit(ind),tau_uncert,'k','LineWidth',2);
%     
% end
% 
% % legend({'Ceres','Linear fit','Ice layer','Homogeneous'},'FontSize',fntsize_sm);
% 
% PrintWhite(plot_relax_times,[fig_folder 'Fig_layer_relaxation_times.jpg']);

%%
% 
% n = 18;
% ind = n + 1;
% 
% figure; hold on;
% 
% cond_order = ((log10(sdl_limb_0(ind)) - log10(sdl_limb(ind,:))) < tol_order);
% first_false_ind = find(~cond_order,1);
% cond_order(first_false_ind:end) = 0;
% 
% [p,S] = polyfit([0 t(find(cond_order)+1)],[log(sdl_limb_0(ind,:)) ...
%     log(sdl_limb(ind,cond_order))],1);
% 
% [log_sdl_fit,log_delta_sdl_fit] = polyval(p,t,S);
% 
% sdl_fit = exp(log_sdl_fit);
% sdl_fit_max = exp(log_sdl_fit+log_delta_sdl_fit);
% sdl_fit_min = exp(log_sdl_fit-log_delta_sdl_fit);
% 
% plot(t,sdl_fit,'ok','LineWidth',3);
% plot(t,sdl_fit_max,'.k','LineWidth',1);
% plot(t,sdl_fit_min,'.k','LineWidth',1);
% 
% plot(0,(sdl_limb_0(ind,:)),'or');
% plot(t(find(cond_order)+1),(sdl_limb(ind,cond_order)),'or');
% 


