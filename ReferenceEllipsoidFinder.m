load sigma_a.mat
% load da.mat
close all;
figure



pcolor(y_surf,z_surf,log10(abs(da*1e5))); shading interp; colorbar;


figure
pcolor(y_surf,z_surf,log10(sigma_a_rad)); shading interp; colorbar;

figure; hold on;
pcolor(y,z,log10((abs(da*1e5)./(sigma_a_rad)))); shading interp;
contour(y,z,log10((abs(da*1e5)./(3*sigma_a_rad))),log10(1*[1 1]),'Color','k','LineWidth',2);
colorbar

hold on

a_ell=292800*1.1;
c_ell=265000*1.1;
x_ell=a_ell*cos(-pi/2:pi/100:pi/2);
z_ell=c_ell*sin(-pi/2:pi/100:pi/2);
plot(x_ell,z_ell,'-b','LineWidth',2);
h_ref_ell=plot(-x_ell,z_ell,'-b','LineWidth',2);

legend(h_ref_ell,{'Reference ellipsoid'},'FontSize',6);