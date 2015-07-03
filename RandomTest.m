ccc;

%% Generate random field with given spectral slope
BTA = -3;
[lmcosi,~,bto,~,el]=plm2rnd(50,BTA,1);

%% Convert to spatial domain
[r,lon,lat] = plm2xyz(lmcosi,0.5);
[lon,lat] = meshgrid(lon,lat);

%% Plot Map
AGUaxes;
pcolorm(lat,lon,r);

%% Compute spectral density
[sdl,l,bta] = plm2spec(lmcosi,3);

%% Plot Spectral density
figure; hold on;
plot(l,sdl,'-ok','LineWidth',3)
set(gca,'FontSize',20);
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Spetral density','FontSize',20);
xlabel('Degree','FontSize',20)
box on

%% Plot 3-D figure
% MeanRadius = 200;
% [x,y,z] = sph2cart(lon/180*pi,lat/180*pi,MeanRadius+r);
% figure; hold on;
% surf(x,y,z,MeanRadius+r); 
% StandardLight

%% Get profile

R=500000;

r_prof = r(:,1); % north-south profile
dist_prof=linspace(0,pi*R,numel(r_prof));
halfcircum = pi*R;
inter = halfcircum/numel(r_prof);
Fs=1/inter;

Y = fft(r_prof);
freq = 0:Fs/length(r_prof):Fs/2;
wavelength=2*pi./freq;
% degree=2*pi*R./wavelength;
k=2*pi*freq;
degree=(sqrt(1+4*k.*k*R*R)-1)/2;

Y = Y(1:length(r_prof)/2+1);
sdl_prof=abs(Y);

% figure; hold on;
% set(gca,'FontSize',20);
% set(gca,'YScale','log');
% set(gca,'XScale','log');
plot(degree,sdl_prof,'-ob','LineWidth',3)

ylim([1e-10 1e10]);
xlim([1e-1 1e2]);


















