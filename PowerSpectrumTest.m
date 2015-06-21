ccc

topo_file1='SH_Vesta20140513shape_64.txt';
topo_file2='SH_VestaShapeDLR64_pixel.txt';


% in=fopen(topo_file,'r');
% fgetl(in)
% fclose(in);

lmcosi1=ReadBalminoSH2(topo_file1);
lmcosi2=ReadBalminoSH2(topo_file2);


%%
% [sdl,l]=plm2spec(lmcosi);

% maxdeg=lmcosi(end,1);
% power=ones(1,maxdeg);
% progressbar(0)
% for n=1:maxdeg     
%     is=n*(n+1)/2+1;
%     ifin=(n+1)*(n+2)/2;   
%     power(n)=sum((lmcosi(is:ifin,3).^2)+(lmcosi(is:ifin,4).^2))/(2*n+1);  
%     progressbar(n/maxdeg);
% end
% progressbar(1);
lmcosi1(:,3:4)=lmcosi1(:,3:4)*1000;

[l1,sdl1]=PowerSpectrum(lmcosi1);
[l2,sdl2]=PowerSpectrum(lmcosi2);


figure; hold on;

plot(1./l1,sdl1,'-r','LineWidth',1);
plot(1./l2,sdl2,'-b','LineWidth',1);
xlabel('Wavelength [km]','FontSize',12);
ylabel('Topography power [km^{3}]','FontSize',12);
set(gca,'YScale','log');
set(gca,'XScale','log');
set(gca,'XDir','reverse');
legend({'SPC','SPG'},'FontSize',12);
box on

%%
figure; hold on;

plot(1./l1,sdl1./sdl2,'-k','LineWidth',1);

xlabel('Wavelength [km]','FontSize',12);
ylabel('Topography power [km^{3}]','FontSize',12);
set(gca,'YScale','log');
set(gca,'XScale','log');
set(gca,'XDir','reverse');
legend({'SPC/SPG'},'FontSize',12);
box on
ylim([0.1/2 10])



%% Isotropic ratio

q1=IsotropicRatio(lmcosi1,lmcosi1);
q2=IsotropicRatio(lmcosi2,lmcosi2);

figure; hold on;
set(gca,'FontSize',12);

plot(1./l1,q1,'r-');
plot(1./l2,q2,'b-');

xlabel('Wavelength [km]','FontSize',12);
ylabel('Isotropuc ratio []','FontSize',12);
set(gca,'YScale','log');
set(gca,'XScale','log');
set(gca,'XDir','reverse');
legend({'SPC','SPG'},'FontSize',12);
box on


%% Correlation

Corr=SphericalHarmonicCorrelation(lmcosi1,lmcosi2);

figure; hold on;
set(gca,'FontSize',12);

plot(1./l1,Corr,'k-');

xlabel('Wavelength [km]','FontSize',12);
ylabel('Correlation []','FontSize',12);
set(gca,'XScale','log');
set(gca,'XDir','reverse');
box on



lmcosi2(mod(lmcosi2(:,1)+lmcosi2(:,2),2)==1,3:4)=...
    lmcosi2(mod(lmcosi2(:,1)+lmcosi2(:,2),2)==1,3:4).*-1;





% figure
% plot(power,'-k','LineWidth',2);
% set(gca,'YScale','log')

% %%
% [r,lon,lat]=plm2xyz(lmcosi,1);
% [lon,lat]=meshgrid(lon/180*pi,lat/180*pi);
%     
% AGUaxes;
% pcolorm(lat,lon,flipud(r)); shading interp
% 
% cbar=colorbar('FontSize',12);
% ylabel(cbar,'Radius [km]','FontSize',12)
% 
% %%
% [x,y,z]=sph2cart(lon,lat,flipud(r));
% 
% figure; hold on;
% 
% surf(x,y,z,flipud(r)); 
% 
% StandardLight