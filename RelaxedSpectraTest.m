ccc;

%% Plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';
%% reading data
% filename_quad = {  '../CeresFE/FE/outer_vertices_0.txt',...
%     '../CeresFE/FE/outer_vertices_1.txt',...
%      '../CeresFE/FE/outer_vertices_2.txt'};

folder_path = '../CeresFE/FE/OutputHomo/';


t = [0:8 10:19]*1e3;

filename_quad = getAllFiles(folder_path);

L = 18;

fig_spec=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('Degree','FontSize',fntsize,'interpreter','latex');
ylabel('Power [$\textrm{km}^{2}$]','FontSize',fntsize,'interpreter','latex');

ylim([1e0 1e8]);
xlim([1 L+10]);


ccj = cool(numel(filename_quad)-1);

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
plot(l_limb_0(1:2:end),sdl_limb_0(1:2:end)...
    ,'-o','MarkerSize',2,'Color','k');
pause(0.5);
drawnow;

for i = 2:numel(filename_quad)
    lmcosi_limb = quad2plm(filename_quad{i},L);
    [sdl_limb(:,i-1),l_limb] = plm2spec(lmcosi_limb);
    
    d = (log10(sdl_limb_0) - log10(sdl_limb(:,i-1)));
    d(3:2:end)
    
    tau(:,i-1) = t(i)./(log(sdl_limb_0) - log(sdl_limb(:,i-1)));
    
    figure(fig_spec);
    pl_relax = plot(l_limb(1:2:end),sdl_limb(1:2:end,i-1)...
        ,'-o','MarkerSize',2,'Color',ccj(i-1,:));
    
    drawnow;
end

PrintWhite(fig_spec,[fig_folder 'Fig_relaxed_spectrum.jpg']);

%% Relxation time
figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'YScale','log');

ccj = cool(size(tau,2));

for i = 1:size(tau,2)
    
    plot(2:2:L,tau(3:2:end,i),'o-k','LineWidth',3,'MarkerSize',5,...
        'Color',ccj(i,:));
    
end

log10(tau(3:2:end,:))

xlabel('Degree','FontSize',fntsize,'interpreter','latex');
ylabel('Characteristic relaxation time [years]','FontSize',fntsize,...
    'interpreter','latex');


% figure;
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on;grid on; box on;

for ind = 3:2:size(sdl_limb,1)
    
    l_fit = l_limb(ind);
    [p,S] = polyfit(t,[log(sdl_limb_0(ind,:)) log(sdl_limb(ind,:))],1);
    tau_fit(ind) = -1/p(1);
    
    cf = fit(t',[log(sdl_limb_0(ind,:)) log(sdl_limb(ind,:))]','poly1');
    
    cf_coeff = coeffvalues(cf);
    cf_confint = confint(cf);
    a = cf_coeff(1);
    b = cf_coeff(2);
    a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
    b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
    
    tau_uncert = sqrt((a_uncert^2)/(a.^4));
    
    plot(ind-1,tau_fit(ind),'o','MarkerSize',6,...
        'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(ind-1,tau_fit(ind),tau_uncert,'k','LineWidth',2);
    
end







