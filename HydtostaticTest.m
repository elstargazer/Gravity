ccc


e=[0.01000 0.02000 0.03000 0.04000 0.05000 0.06000 0.07000 0.08000 0.09000 0.10000 0.11000 0.12000 0.13000 0.14000 0.15000 0.16000 0.17000 0.18000 0.19000 0.20000 0.21000 0.22000 0.23000 0.24000];
E1=[0.13959 0.19761 0.24158 0.27886 0.31147 0.34050 0.36770 0.39250 0.41592 0.43790 0.45864 0.47895 0.49789 0.51604 0.53399 0.55103 0.56739 0.58335 0.59873 0.61401 0.62897 0.64315 0.65732 0.67128];
E2=[0.14390 0.20330 0.24860 0.28670 0.32010 0.35010 0.37770 0.40320 0.42710 0.44960 0.47090 0.49130 0.51070 0.52930 0.54730 0.56460 0.58130 0.59750 0.61320 0.62854 0.64350 0.65800 0.67224 0.68620];
E2RD=[0.14383 0.20288 0.24782 0.28540 0.31824 0.34769 0.37453 0.39931 0.42238 0.44402 0.46441 0.48372 0.50208 0.51957 0.53630 0.55232 0.56771 0.58250 0.59674 0.61047 0.62373 0.63654 0.64894 0.66093];

G=6.67e-11;
rho2=3900;
T=sqrt(2*pi./e/G/rho2);

f2=1-(1-E2.^2).^(0.5);
f1=1-(1-E1.^2).^(0.5);
f2RD=1-(1-E2RD.^2).^(0.5);


f2TOF=(-3.65657 + 3.02014*(T/3600).^2)./((T/3600).^4);


figure; hold on;
set(gca,'FontSize',20);

plot(T/3600,f1./f1,'o--b')
plot(T/3600,f2./f1,'o-r');
plot(T/3600,f2RD./f1,'o-g');
plot(T/3600,f2TOF./f1,'o-m');


% figure; hold on;
% 
% plot(e,E1,'-b')
% plot(e,E2RD,'--r');

legend({'core flattening exact','Outer shape flattening exact','Outer shape flattening Radau-Darwin ','Outer shape Dermott 1979'},'FontSize',20);
xlabel('Period [hours]','FontSize',20);
ylabel('Flattening','FontSize',20);

xlim([min(T/3600) max(T/3600)]);

