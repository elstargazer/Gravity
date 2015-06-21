MaxDegreeTopo=20;
load VestaHASTALAVESTAshape_sh720.mat
lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

MaxTopoPower=20;
MaxDegreeTopo=20;
MaxDegreeGrav=20;

ro=3.458942584995802e+03;
G=6.67384e-11;
mu=17.2883080293e9;
Rref=265000;

Resolution=1;

[ri_shape,~,~,~]=plm2xyz(lmcosi_shape,Resolution);

lmcosi_gt=TopoSH2GravitySH(ri_shape,mu,ro,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);



GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_g,Rref,mu_g,mu_std]=ReadGRAILGravityModel(GravityFileName);

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20F/JGV20F02.SHA';
[lmcosi_f,Rref,mu_f,mu_std]=ReadGRAILGravityModel(GravityFileName);

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20E/JGV20E02.SHA';
[lmcosi_e,Rref,mu_e,mu_std]=ReadGRAILGravityModel(GravityFileName);

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA10/JGV10C02.SHA';
[lmcosi_c,Rref,mu_c,mu_std]=ReadGRAILGravityModel(GravityFileName);

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA6/JGV06B02.SHA';
[lmcosi_b,Rref,mu_b,mu_std]=ReadGRAILGravityModel(GravityFileName);

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA6/JGV06A01.SHA';
[lmcosi_a,Rref,mu_a,mu_std]=ReadGRAILGravityModel(GravityFileName);



% [sdl_a,l_a,~,~,~,~]=plm2spec(lmcosi_a,3);
% [sdl_a_error,l_a_error,~,~,~,~]=plm2spec([lmcosi_a(:,1:2) lmcosi_a(:,5:6)],3);

[sdl_b,l_b,~,~,~,~]=plm2spec(lmcosi_b);
[sdl_b_error,l_b_error,~,~,~,~]=plm2spec([lmcosi_b(:,1:2) lmcosi_b(:,5:6)]);

[sdl_c,l_c,~,~,~,~]=plm2spec(lmcosi_c,3);
[sdl_c_error,l_c_error,~,~,~,~]=plm2spec([lmcosi_c(:,1:2) lmcosi_c(:,5:6)],3);

[sdl_e,l_e,~,~,~,~]=plm2spec(lmcosi_e,3);
[sdl_e_error,l_e_error,~,~,~,~]=plm2spec([lmcosi_e(:,1:2) lmcosi_e(:,5:6)],3);

[sdl_f,l_f,~,~,~,~]=plm2spec(lmcosi_f,3);
[sdl_f_error,l_f_error,~,~,~,~]=plm2spec([lmcosi_f(:,1:2) lmcosi_f(:,5:6)],3);

[sdl_g,l_g,~,~,~,~]=plm2spec(lmcosi_g,3);
[sdl_g_error,l_g_error,~,~,~,~]=plm2spec([lmcosi_g(:,1:2) lmcosi_g(:,5:6)],3);

[figure1,axes1]=WhiteFigureLog([],'Coefficient Magnitude','Degree');

% h=plot(l_a,sdl_a,l_a_error,sdl_a_error,...
%     l_c,sdl_c,l_c_error,sdl_c_error,...
%     l_e,sdl_e,l_e_error,sdl_e_error,...
%     l_f,sdl_f,l_f_error,sdl_f_error,...
%     l_g,sdl_g,l_g_error,sdl_g_error);

h1 =plot(l_b,sdl_b,'or-','LineWidth',2,'MarkerSize',3);
h2 =plot(l_c,sdl_c,'og-','LineWidth',2,'MarkerSize',3);
h3 =plot(l_e,sdl_e,'ob-','LineWidth',2,'MarkerSize',3);
% h4=plot(l_f,sdl_f,'om-',l_f_error,sdl_f_error,'om-');
h5= plot(l_g,sdl_g,'om-','LineWidth',2,'MarkerSize',3);

semilogy(l_bouguer,sdl_bouguer,'o-','Color',[0.2 0.4 0.5],'LineWidth',2,'MarkerSize',3);
semilogy(l_bouguer_min,sdl_bouguer_min,'o-','Color',[0.7 0.5 0.3],'LineWidth',2,'MarkerSize',3);



kaula_law=0.011./(l_g.^2).*sqrt(2*l_g+1);
kaula_law_plot=plot(l_g,kaula_law,'ok-','LineWidth',2,'MarkerSize',3);


[sdl_gt,l_gt,~,~,~,~]=plm2spec(lmcosi_gt,3);
homogenen_plot=plot(l_gt,sdl_gt,'oc-','LineWidth',2,'MarkerSize',3);

h1e=plot(l_b_error,sdl_b_error,'or--','LineWidth',2,'MarkerSize',3);
h2e=plot(l_c_error,sdl_c_error,'og--','LineWidth',2,'MarkerSize',3);
h3e=plot(l_e_error,sdl_e_error,'ob--','LineWidth',2,'MarkerSize',3);
h5e=plot(l_g_error,sdl_g_error,'om--','LineWidth',2,'MarkerSize',3);



% legend({'JGV06A01','JGV06A01 error',...
%     'JGV10C02','JGV10C02 error',...
%     'JGV20E02','JGV20E02 error',...
%     'JGV20F02','JGV20F02 error',...
%     'JGV20G02','JGV20G02 error'},'TextColor','k')


set(axes1,'XTick',2:2:20);
set(axes1,'YTick',10.^(-8:-1));

axis([2 20 1e-9 1e-1]);


% set(gca, 'Position',[0 0 1 1])
set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

legend({'Survey (vesta6b)','HAMO1 (vesta10c)','LAMO (vesta20f)','HAMO2 (vesta20g)','Bouguer','Min Bouguer','Kaula law','Homogeneous density'},'TextColor','k','FontSize',8,'Location','SouthEast')

print(figure1, '-dpsc', 'GravitySpectrum.eps');

