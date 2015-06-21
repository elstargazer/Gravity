ccc
load Shape1800

MaxDegreeTopo=150;
lmcosi_dt=TruncateGravityModel(lmcosi_d,MaxDegreeTopo,1);
lmcosi_gt=TruncateGravityModel(lmcosi_g,MaxDegreeTopo,1);


%% Finding initial correlation
Corr=SphericalHarmonicCorrelation(lmcosi_dt,lmcosi_gt);


%% Plotting initial correlation
figure; hold on;
set(gca,'FontSize',20);
plot(Corr,'b','LineWidth',3);



%% Searching minimum

% alp0=0.001;
% bta0=0.001;
% gam0=0.001;


options=optimset('TolX',.02/3600);


alp0=300*(rand(1,8)-0.5)/3600;
bta0=300*(rand(1,8)-0.5)/3600;
gam0=300*(rand(1,8)-0.5)/3600;

x0=[alp0; bta0; gam0];

x=x0;

parfor i=1:8

[x(:,i),fval(i),~] = fminsearch(@(ang) -SphCorrRot(lmcosi_dt,lmcosi_gt,ang),x0(:,i),options);

end

ind=find(fval==min(fval),1);

%% Computing final correlation

[lmcosi_dtrot,~,~]=plm2rot(lmcosi_dt,x(1,ind),x(2,ind),x(3,ind));
Corr2=SphericalHarmonicCorrelation(lmcosi_dtrot,lmcosi_gt);

disp(x*3600);

%% Plotting  final correlation
plot(Corr2,'r','LineWidth',3);