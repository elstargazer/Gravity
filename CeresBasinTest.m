ccc

%% Input parameters
L = 5;
G = 6.67e-11;

GM=62.68e9; % PCK version 0.2
T=9.073859324514187; % DLR

a=484170;
b=481410;
c=447800;

M=GM/G;
R=(a.*a.*c).^(1/3);
V=1e9*4/3*pi*R.^3;

router = R;

aref = 482000+350000;
cref = 446000+350000;

Rref_grav = 482000;
MaxDegreeTopo = 50;
MaxDegreeGrav = 20;
MaxTopoPower = 5;

rcore=300000:1500:470000;
rhocore=2000:500:6000;

h = 25000;
d = 100000;

[rhocorei,rcorei]=meshgrid(rhocore,rcore);
rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));

figure; hold on;
pcolor(rcorei,rhocorei,rhoouteri); shading interp
caxis([0 4000]);

i = 200;
plot(rcorei(i), rhocorei(i),'wo');

%% make a mesh
step = 1;
fi = (-90:step:90)/180*pi;
lambda = (-180:step:180)/180*pi;
[lambda, fi] = meshgrid(lambda, fi);

% basin loction
fic=90/180*pi;
lambdac=0/180*pi;

%% compute hydrostatic equilibrium shape

[fhi,fval]=HydrostaticStateExact2l(router,rcorei(i),T,rhoouteri(i),rhocorei(i),0.1,0.1);

fp1 = fhi(1);
fp2 = fhi(2);

%% setup moho surface

[acore,bcore,ccore]=fr2abc(rcorei(i),fp1,0);

delta_rho = rhocorei(i) - rhoouteri(i);

[xcore, ycore, zcore] = TriEllRadVec(fi,lambda,acore,bcore,ccore,'xyz');
rcore = sqrt(xcore.*xcore+ycore.*ycore+zcore.*zcore);

dist = distance(fic,lambdac,fi,lambda,...
    [sqrt(acore*bcore) Eccentricity(sqrt(acore*bcore),ccore)],...
    'radians');

uplift = h.*exp(-(dist/d).^2);

figure; hold on; 
set(gca,'FontSize',20);
box on;
plot(dist(:,1)/1000,uplift(:,1)/1000,'-k','LineWidth',3);
xlabel('Distance [km]','FontSize',20);
ylabel('Uplift [km]','FontSize',20);

Vcore = 4/3*pi*acore*bcore*ccore;
mucore = Vcore*(rhocorei(i) - rhoouteri(i))*G;

lmcosi_g_core=TopoSH2GravitySH(rcore, mucore,rhocorei(i),Rref_grav,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

rcore = rcore + uplift;

[xcore, ycore, zcore] = sph2cart(lambda,fi,rcore);

Vcore = Mesh2Volume(xcore,ycore,zcore);
mucore_up = Vcore*(rhocorei(i) - rhoouteri(i))*G;

lmcosi_g_core_up=TopoSH2GravitySH(rcore, mucore_up,rhocorei(i),Rref_grav,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

%% setup outer surface

[x,y,z] = TriEllRadVec(fi,lambda,a,b,c,'xyz');
r = sqrt(x.*x + y.*y + z.*z);

Vouter = Mesh2Volume(x,y,z);

muouter = Vouter*rhoouteri(i)*G;

lmcosi_g_outer=TopoSH2GravitySH(r,muouter,rhoouteri(i),Rref_grav,...
    MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

%% draw structure

figure; hold on;

surf(xcore,ycore,zcore,rcore);
surf(x,y,z,r);

colormap jet
axis equal; shading interp; alpha(0.5);

%% compute gravity coefficients

lmcosi_g_total = lmcosi_g_outer;
lmcosi_g_total_up = lmcosi_g_outer;
mutotal_up = mucore_up + muouter;
mutotal = mucore + muouter;

lmcosi_g_total(:,3:4) = (lmcosi_g_outer(:,3:4) * muouter + ...
    lmcosi_g_core(:,3:4)*mucore)/mutotal;

lmcosi_g_total_up(:,3:4) = (lmcosi_g_outer(:,3:4) * muouter + ...
    lmcosi_g_core_up(:,3:4)*mucore_up)/mutotal_up;

%% plot gravity 

[xsurf,ysurf,zsurf] = TriEllRadVec(fi,lambda,aref,aref,cref,'xyz');

lmcosi_g_total = TruncateGravityModel(lmcosi_g_total,L,1);
lmcosi_g_total_up = TruncateGravityModel(lmcosi_g_total_up,L,1);

[gx,gy,gz] = GravityAcceleration(mutotal,Rref_grav,lmcosi_g_total,...
    xsurf,ysurf,zsurf);
[gx_up,gy_up,gz_up] = GravityAcceleration(mutotal,Rref_grav,lmcosi_g_total_up,...
    xsurf,ysurf,zsurf);
[g_up,~,~]=GravityComponents(gx,gy,gz,...
    xsurf,ysurf,zsurf,aref,cref);
[g_up_up,~,~]=GravityComponents(gx_up,gy_up,gz_up,...
    xsurf,ysurf,zsurf,aref,cref);

%% gravity profile

anom = (g_up_up - g_up)*1e5;

AGUaxes;
pcolorm(fi*180/pi,lambda*180/pi,anom); shading interp;
cbar = colorbar('FontSize',20);
ylabel(cbar,'Gravity anomaly [mGal]','FontSize',20);

figure; hold on;
plot(fi(:,1),anom(:,1),'-k','LineWidth',3)

figure; hold on;
plot(fi(:,1),g_up(:,1)*1e5,'-k','LineWidth',3)
