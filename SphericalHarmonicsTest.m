% ccc

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';

N_trunc=15;

[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);
lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);

% lmcosi_obs(1:end-1,3:4)=0;

r=300000;
lambda=0:1:360;
fi=-90:1:90;

[lambda,fi]=meshgrid(lambda,fi);
% fi=180*rand(1,100)-90;

[x,y,z]=sph2cart(lambda/180*pi,fi/180*pi,r);

U1=GravityPotential2(mu,Rref,lmcosi_obs,x,y,z);

lmcosi_pot=plm2pot(AddZeroHarm(lmcosi_obs,1),r,mu,Rref,1,'nothing');

[U2,~,~,~]=plm2xyz(lmcosi_pot,fi(:),lambda(:));
U2=reshape(U2,size(x));

AGUaxes;
set(gca,'FontSize',20);

d=(U1-U2)./U1;

pcolorm(fi/180*pi,lambda/180*pi,d); colorbar;

% U2=U2-mu/r
% Um=1.6134523131878e-6;
% [U1; U2; Um]

% plot((U1-U2)./U1)

