ccc

EarthFileName='/Users/antonermakov/Earth/srtmp2160.ret_shape';
MoonFileName='/Users/antonermakov/GRAIL/Topography/SH/LRO_LTM05_2050_SHA.TAB.txt';
VenusFileName='/Users/antonermakov/Venus/VenusTopo719.shape';
MarsFileName='/Users/antonermakov/Mars/MarsTopo719.shape';

% VestaMatFileName='VestaHASTALAVESTAshape_sh1000.mat';
VestaMatFileName='~/Dawn/SH/AnalysisSynthesis/VestaTest/SH_VestaHASTALAVESTAshape_6min';

%% Vesta topo sh

% load(VestaMatFileName);lmcosi=lmcosi_shape;
lmcosi_Vesta=ReadBalminoSH2(VestaMatFileName);

%% Planets topo sh
lmcosi_Moon=load(MoonFileName);
lmcosi_Venus=load(VenusFileName);
lmcosi_Mars=load(MarsFileName);
lmcosi_Earth=load(EarthFileName);

lmcosi_Ceres(:,3:4)=lmcosi_Ceres(:,3:4)*1000;

% lmcosi_Moon(:,3)=lmcosi_Moon(:,3)/lmcosi_Moon(1,3);
% lmcosi_Moon(:,4)=lmcosi_Moon(:,4)/lmcosi_Moon(1,3);
% 
% lmcosi_Venus(:,3)=lmcosi_Venus(:,3)/lmcosi_Venus(1,3);
% lmcosi_Venus(:,4)=lmcosi_Venus(:,4)/lmcosi_Venus(1,3);
% 
% lmcosi_Mars(:,3)=lmcosi_Mars(:,3)/lmcosi_Mars(1,3);
% lmcosi_Mars(:,4)=lmcosi_Mars(:,4)/lmcosi_Mars(1,3);
% 
% lmcosi_Earth(:,3)=lmcosi_Earth(:,3)/lmcosi_Earth(1,3);
% lmcosi_Earth(:,4)=lmcosi_Earth(:,4)/lmcosi_Earth(1,3);
% 
%  lmcosi(:,3)=lmcosi(:,3)/lmcosi(1,3);
%  lmcosi(:,4)=lmcosi(:,4)/lmcosi(1,3);


% [sdl_Moon,l_Moon,bta_Moon,lfit_Moon,logy_Moon,logpm_Moon]=plm2spec(lmcosi_Moon);
% [sdl_Venus,l_Venus,bta_Venus,lfit_Venus,logy_Venus,logpm_Venus]=plm2spec(lmcosi_Venus);
% [sdl_Mars,l_Mars,bta_Mars,lfit_Mars,logy_Mars,logpm_Mars]=plm2spec(lmcosi_Mars);
% [sdl_Vesta,l_Vesta,bta_Vesta,lfit_Vesta,logy_Vesta,logpm_Vesta]=plm2spec(lmcosi);

[k_Moon,sdl_Moon]=PowerSpectrum(lmcosi_Moon);
[k_Mars,sdl_Mars]=PowerSpectrum(lmcosi_Mars);
[k_Venus,sdl_Venus]=PowerSpectrum(lmcosi_Venus);
[k_Earth,sdl_Earth]=PowerSpectrum(lmcosi_Earth);
[k_Vesta,sdl_Vesta]=PowerSpectrum(lmcosi_Vesta);
[k_Ceres,sdl_Ceres]=PowerSpectrum(lmcosi_Ceres);

figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],...
    'YTick',10.^(0:2:14),...
    'YScale','log',...
    'YMinorTick','on',...
    'YMinorGrid','on',...
    'YColor',[0 0 0],...
    'XTick',[1e-05 0.0001 0.001 0.01 0.1 1],...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'XColor',[0 0 0],...
    'FontSize',12,...
    'FontName','Helvetica',...
    'Color',[1 1 1]);
box(axes1,'on');
hold(axes1,'all');
hold on;
% Create semilogy
% loglog(l_Moon,sdl_Moon,'r','MarkerFaceColor','r','MarkerSize',1,'Marker','o',...
%      'LineWidth',1);
% loglog(l_Venus,sdl_Venus,'g','MarkerFaceColor','g','MarkerSize',1,'Marker','o',...
%  'LineWidth',1);
% loglog(l_Mars,sdl_Mars,'b','MarkerFaceColor','b','MarkerSize',1,'Marker','o',...
%     'LineWidth',1);
% loglog(l_Vesta,sdl_Vesta,'y','MarkerFaceColor','b','MarkerSize',1,'Marker','o',...
%     'LineWidth',1);

loglog(k_Moon,sdl_Moon,'r','MarkerFaceColor','r','MarkerSize',1,'Marker','o',...
     'LineWidth',1);
loglog(k_Venus,sdl_Venus,'g','MarkerFaceColor','g','MarkerSize',1,'Marker','o',...
 'LineWidth',1);
loglog(k_Mars,sdl_Mars,'b','MarkerFaceColor','b','MarkerSize',1,'Marker','o',...
    'LineWidth',1);
loglog(k_Earth,sdl_Earth,'c','MarkerFaceColor','b','MarkerSize',1,'Marker','o',...
    'LineWidth',1);
loglog(k_Vesta,sdl_Vesta,'m','MarkerFaceColor','b','MarkerSize',1,'Marker','o',...
    'LineWidth',1);
loglog(k_Ceres,sdl_Ceres,'k','MarkerFaceColor','k','LineWidth',2,'MarkerSize',1,'Marker','o',...
    'LineWidth',1);



% legend({['Earth' num2str(10.^p_Earth(2)) '* k^' num2str(p_Earth(1)],...
%     ['Venus' num2str(10.^p_Venus(2)) '* k^' num2str(p_Venus(1)],...
%     ['Mars' num2str(10.^p_Mars(2)) '* k^' num2str(p_Mars(1)],...
%     ['Moon' num2str(10.^p_Moon(2)) '* k^' num2str(p_Moon(1)]...,
%     ['Vesta' num2str(10.^p_Vesta(2)) '* k^' num2str(p_Vesta(1)]}...
%     ,'TextColor','w');


% loglog(lfit_Moon,logpm_Moon,'--r','LineWidth',1);
% loglog(lfit_Venus,logpm_Venus,'--g','LineWidth',1);
% loglog(lfit_Mars,logpm_Mars,'--b','LineWidth',1);
% loglog(lfit_Vesta,logpm_Vesta,'--y','LineWidth',1);

% title({['Moon spectral slope = ' num2str(bta_Moon)],...
%     ['Venus spectral slope = ' num2str(bta_Venus)],...
%     ['Mars spectral slope = ' num2str(bta_Mars)],...
%     ['Vesta spectral slope = ' num2str(bta_Venus)]},...
%     'fontsize',25,'Color','w');


grid on;

ylabel('Total power [km^3]');
xlabel('k [cycles/km]');

% title('Energy spectral density of topography','Color','w');

p_Earth = polyfit(log10(k_Earth),log10(sdl_Earth),1);
p_Moon = polyfit(log10(k_Moon),log10(sdl_Moon),1);
p_Venus = polyfit(log10(k_Venus),log10(sdl_Venus),1);
p_Mars = polyfit(log10(k_Mars),log10(sdl_Mars),1);
p_Vesta = polyfit(log10(k_Vesta),log10(sdl_Vesta),1);
p_Ceres = polyfit(log10(k_Ceres),log10(sdl_Ceres),1);

% loglog(k_Earth,(10.^p_Earth(2))*(k_Earth.^p_Earth(1)),'-c','markersize',0.1);
% loglog(k_Venus,(10.^p_Venus(2))*(k_Venus.^p_Venus(1)),'-g','markersize',0.1);
% loglog(k_Mars,(10.^p_Mars(2))*(k_Mars.^p_Mars(1)),'-b','markersize',0.1);
% loglog(k_Moon,(10.^p_Moon(2))*(k_Moon.^p_Moon(1)),'-r','markersize',0.1);
% loglog(k_Vesta,(10.^p_Vesta(2))*(k_Vesta.^p_Vesta(1)),'-y','markersize',0.1);


% legend({['Moon \approx ' num2str(10.^p_Moon(2),'%3.2e')  '\cdotk^{' num2str(p_Moon(1),'%3.2f') '}' ', D = ' num2str((5+p_Moon(1))/2,'%3.2f')],...
%     ['Venus \approx ' num2str(10.^p_Venus(2),'%3.2e') '\cdotk^{' num2str(p_Venus(1),'%3.2f') '}' ', D = ' num2str((5+p_Venus(1))/2,'%3.2f')],...
%     ['Mars \approx  ' num2str(10.^p_Mars(2),'%3.2e') '\cdotk^{' num2str(p_Mars(1),'%3.2f') '}' ', D = ' num2str((5+p_Mars(1))/2,'%3.2f')],...
%     ['Earth \approx ' num2str(10.^p_Earth(2),'%3.2e') '\cdotk^{' num2str(p_Earth(1),'%3.2f') '}' ', D = ' num2str((5+p_Earth(1))/2,'%3.2f')]...,
%     ['Vesta \approx ' num2str(10.^p_Vesta(2),'%3.2e') '\cdotk^{' num2str(p_Vesta(1),'%3.2f') '}' ', D = ' num2str((5+p_Vesta(1))/2,'%3.2f')] }...
%     ,'TextColor','k','fontsize',12);

legend({'Moon','Venus','Mars','Earth','Vesta','Ceres'},'FontSize',12,'Location','SouthWest');

set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

print(figure1, '-dpsc', 'TopographySpectrum.eps');



