filename_g='SH_Vesta20130522shape_6m_mean';
filename_d='SH_VestaShapeDLR_6m_mean';

lmclmcosi_d=ReadBalminoSH2(filename_d);

[sdl_g,l_g]=plm2spec(lmcosi_g,3);
[sdl_d,l_d]=plm2spec(lmcosi_d,3);




%% Rotate

[lmcosi_g,spec1,spec2]=plm2rot(lmcosi_g,0,0,0.0556,'dlmb',1);

[ri_g,lambda,fi]=plm2xyz(lmcosi_g,.1); ri_g=ri_g*1e3;


[ri_d,~,~]=plm2xyz(lmcosi_d,.1);

[lambdai,fii]=meshgrid(lambda,fi);

% figure
% pcolor(ri_d);shading interp
% 
% 
% figure
% pcolor(ri_g);shading interp
% 
% 
% figure;
% surf(((ri_g-ri_d)));shading interp; view(0,90)

% caxis([-400 400])

%% Plot Power
figure; hold on;
set(gca,'FontSize',20);
set(gca,'YScale','log');

semilogy(l_g,sdl_g*1e3,'b')
semilogy(l_d,sdl_d,'r')

xlabel('Degree','FontSize',20);
ylabel('Sum of Coeffs^2','FontSize',20);
legend({'Gaskell','DLR'},'FontSize',20);

%% Plot Differnce of Power
figure; hold on;
semilogy(l_g,(sdl_g*1e3-sdl_d)./sdl_d,'k')
% semilogy(l_g,-(sdl_g*1e3-sdl_d)./sdl_d,'b')
set(gca,'FontSize',20);
% set(gca,'YScale','log');
xlabel('Degree','FontSize',20);
ylabel('Difference of Sum of Coeffs^2','FontSize',20);

%% Plot Power of Difference
lmcosi_diff=lmcosi_g;

lmcosi_diff(:,3)=lmcosi_g(:,3)*1e3-lmcosi_d(:,3);
lmcosi_diff(:,4)=lmcosi_g(:,4)*1e3-lmcosi_d(:,4);


[sdl_diff,~]=plm2spec(lmcosi_diff,3);

figure; hold on;
semilogy(l_g,sdl_diff./sdl_d,'k')
set(gca,'FontSize',20);
% set(gca,'YScale','log');
xlabel('Degree','FontSize',20);
ylabel('Difference of Sum of Coeffs^2','FontSize',20);



%% To XYZ
[x_g,y_g,z_g]=sph2cart(lambdai/180*pi,fii/180*pi,ri_g);
[x_d,y_d,z_d]=sph2cart(lambdai/180*pi,fii/180*pi,ri_d);


%% To geodetic
at=280869.77;
bt=226491.24;
et=Eccentricity(at,bt);

[~,~,H_g]=XYZ2BLH(x_g,y_g,z_g,at,et);
[~,~,H_d]=XYZ2BLH(x_d,y_d,z_d,at,et);

%% Plot 

MapRadialGrid(flipud(H_g))
MapRadialGrid(flipud(H_d))
MapRadialGrid(H_g-H_d)
caxis([-400 400])

%% Correlation
figure
SphericalHarmonicCorrelation(lmcosi_g,lmcosi_d,'r','-');

%% Rotation 
% d=0.052:0.0001:0.06
d=-0.01:0.0001:0.02;

for i=1:numel(d)
    
H_gi=interp1(lambda,H_g(180,:),lambda+d(i));
Con=isnan(H_gi);
H_gi(Con)=0;
H_d(180,Con)=0;
V(i)=sum((H_gi-H_d(180,:)).^2);

end

figure;
plot(d,V)




