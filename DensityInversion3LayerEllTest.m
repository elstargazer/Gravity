ccc

%% Input values

% input gravity and shape model files
GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTAELL/nomfil_VESEH20D3.txt';
filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';

NTess=6;
NTessPlus=-3;
k=3;

if (k<-NTessPlus)
    error('Not enough data points');
end

N_trunc=14;

% mantle
a2=257000;
c2=207000;

% core
a3=117000;
c3=105000;

% reference ellipsoid
aref=303860;
bref=288660;
cref=246460;
RefEll=[aref bref cref];

Re=265000;

% mean density
rho_mean=3457;

% internal structure densities
rho_core=7400;
rho_mantle=3160;
rho_crust=2885;

% gridstep for plotting residuals
FiStep=5;
LambdaStep=5;

% core offset
% x3c=2.2291e+03;
% y3c=9.4132e+03;
% z3c=-29.7072;h = sqrt(a^2-b^2);

mu=17.28830802930000e9;

angle_fi=1.209/180*pi;
angle_theta=-0.558/180*pi;
angle_psi=129.209/180*pi;

% core at the center of the northern ellipsoid
% x3cm=-830;
% y3cm=-200;
% z3cm=-5660;

% make no com-cof offset
x3cm=1.6313e+03;
y3cm=6.8891e+03;
z3cm=-21.7411;

% no core offset
% x3cm=0;
% y3cm=0;
% z3cm=0;

% shape offset
dx=0;
dy=0;
dz=0;

%% Load gravity model

lmcosi=ReadEllHarmModel(GravityFileName);
lmcosi=TruncateEllModel(lmcosi,N_trunc);

%% Load topography model

[TR,x1,y1,z1]=Grid2TriShapeModel(filename,NTess);
[TRo,x1o,y1o,z1o]=Grid2TriShapeModel(filename,NTess-k);
FV=TR.Triangulation;
FVo=TRo.Triangulation;

x1=x1*1000;
y1=y1*1000;
z1=z1*1000;

s=size(FV);
[lambdai,fii,r1]=cart2sph(x1,y1,z1);

% MaxDegreeTopo=7;
% load VestaHASTALAVESTAshape_sh720.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
%
% r1sh=plm2xyz(lmcosi_shape,fii*180/pi,lambdai*180/pi);
% [x1,y1,z1]=sph2cart(lambdai,fii,r1sh);

x1m=(x1(FV(:,1))+x1(FV(:,2))+x1(FV(:,3)))/3;
y1m=(y1(FV(:,1))+y1(FV(:,2))+y1(FV(:,3)))/3;
z1m=(z1(FV(:,1))+z1(FV(:,2))+z1(FV(:,3)))/3;

[lambda1m,fi1m,~]=cart2sph(x1m,y1m,z1m);

%% Create reference surface

TRr=IcosahedronMesh;
TR_2r=SubdivideSphericalMesh(TRr,NTess+NTessPlus);
figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FVr=TR_2r.Triangulation;
xr=aref*TR_2r.X(:,1);
yr=aref*TR_2r.X(:,2);
zr=cref*TR_2r.X(:,3);
xr=(xr(FVr(:,1))+xr(FVr(:,2))+xr(FVr(:,3)))/3;
yr=(yr(FVr(:,1))+yr(FVr(:,2))+yr(FVr(:,3)))/3;
zr=(zr(FVr(:,1))+zr(FVr(:,2))+zr(FVr(:,3)))/3;

clear FVr TR_2r TRr

% amp=5000;
% xr=xr+amp*rand(size(xr));
% yr=yr+amp*rand(size(yr));
% zr=zr+amp*rand(size(zr));

[lambdar,fir,~]=cart2sph(xr,yr,zr);


[xr,yr,zr]=TriEllRadVec(fir,lambdar,aref,bref,cref);

% for fast plotting of united triangles
TRi=IcosahedronMesh;
TRi=SubdivideSphericalMesh(TRi,NTess-k);

xi=aref*TRi.X(:,1);
yi=aref*TRi.X(:,2);
zi=cref*TRi.X(:,3);


%% Plot

figure; hold on; axis equal
tri_fig=trisurf(FV,x1,y1,z1,r1); shading interp
set(tri_fig,'edgecolor','k');
trimesh_fig=trimesh(FV,x1,y1,z1,r1,'edgecolor','k','facecolor','none');

% figure; hold on; axis equal
% tri_fig=trisurf(FV,x1sh,y1sh,z1sh,r1sh); shading interp
% set(tri_fig,'edgecolor','k');
% trimesh_fig=trimesh(FV,x1sh,y1sh,z1sh,r1sh,'edgecolor','k','facecolor','none');

% surf(xr,yr,zr,a);

plot3(xr,yr,zr,'.k','MarkerSize',10);
axis(1.4*[-300000 300000 -300000 300000 -300000 300000])
alpha(0.2);

%% Constructing internal structure

% x2=a2*xu;
% y2=a2*yu;
% z2=c2*zu;
%
% x3=a3*xu;
% y3=a3*yu;
% z3=c3*zu;

[x2,y2,z2]=TriEllRadVec(fii,lambdai,a2,a2,c2);
% [x3,y3,z3]=TwoEllRadVec(lambdai,fii,a3,c3);
[x3,y3,z3]=TriEllShiftRadVec(fii,lambdai,a3,a3,c3,x3cm,y3cm,z3cm);

r2=sqrt(x2.*x2+y2.*y2+z2.*z2);
r3=sqrt(x3.*x3+y3.*y3+z3.*z3);

%
% figure; hold on; axis equal
% tri_fig1=trisurf(FV,x1,y1,z1,r1,'edgecolor','k');
% tri_fig2=trisurf(FV,x2,y2,z2,r2,'edgecolor','k');
% tri_fig3=trisurf(FV,x3,y3,z3,r3);
%
%  shading interp;
% alpha(0.5);

rho1=ones(1,s(1))*rho_crust;
rho2=ones(1,s(1))*rho_mantle-rho1;
rho3=ones(1,s(1))*rho_core-rho2-rho1;

% rotate and shift shape

x1=x1-dx;
y1=y1-dy;
z1=z1-dz;

x2=x2-dx;
y2=y2-dy;
z2=z2-dz;

x3=x3-dx;
y3=y3-dy;
z3=z3-dz;

S=rot(angle_psi,3)*rot(angle_theta,2)*rot(angle_fi,1);

[x1,y1,z1]=RotateByMatrix(x1,y1,z1,S);
[x2,y2,z2]=RotateByMatrix(x2,y2,z2,S);
[x3,y3,z3]=RotateByMatrix(x3,y3,z3,S);
% [xi,yi,zi]=RotateByMatrix(xi,yi,zi,S);
[x1m,y1m,z1m]=RotateByMatrix(x1m,y1m,z1m,S);


%% Computing gravity
% [gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_obs,xr,yr,zr);
% [g_up,g_east,g_north]=GravityComponents(gx,gy,gz,xr,yr,zr,aref,cref);


% generating observed potential
% rho1_i=ones(1,s(1))*rho_crust;
% rho1_i=rho1_i+200*randn(size(rho1_i));
% rho2_i=ones(1,s(1))*rho_mantle-rho1_i;
% rho3_i=ones(1,s(1))*rho_core-rho2_i-rho1_i;
%
%
% U1_o=GravityPotentialTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho1_i);
% U2_o=GravityPotentialTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho2_i);
% U3_o=GravityPotentialTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho3_i);
%
% U_o=U1_o+U2_o+U3_o;

% tic
% [U_o,~,~,~]=plm2xyz(lmcosi_pot,fir*180/pi,lambdar*180/pi);
% t_pot_obs=toc

% U_o=GravityPotential2(mu,Rref,lmcosi_obs,xr,yr,zr);
[U_o,lambda_ell]=GravityPotentialEll4(mu,aref,RefEll,lmcosi,xr,yr,zr);

U_o=U_o/1000;

% gtot=sqrt(gx.*gx+gy.*gy+gz.*gz);

% [ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% atot=sqrt(ax.*ax+ay.*ay+az.*az);
% U_c=GravityPotentialTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% d2=DeltaSqGravity(rho,xi,yi,zi,FV,xr,yr,zr,gx,gy,gz)

%% Numerical Jacobian

% Jacobian linear
% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho,refell);
% a_up0=a_up;

tic
U1_c=GravityPotentialTriDen93(x1',y1',z1',FV,xr,yr,zr,rho1);
U2_c=GravityPotentialTriDen93(x2',y2',z2',FV,xr,yr,zr,rho2);
U3_c=GravityPotentialTriDen93(x3',y3',z3',FV,xr,yr,zr,rho3);
t_pot_calc=toc

U_c=U1_c+U2_c+U3_c;

clear U1_c U2_c U3_c

Nelems=numel(rho1);
Nobs=numel(xr);

% jac=zeros(Nobs,Nelems);

rho1=ones(1,s(1))*rho_crust;
rho2=ones(1,s(1))*rho_mantle-rho1;
rho3=ones(1,s(1))*rho_core-rho2-rho1;

% progressbar(0);
% delta_rho=50;
% J1n=zeros(Nobs,Nelems);
% J2n=zeros(Nobs,Nelems);
%
% for i=1:Nelems   % computing jacobian
%
%     rho1_inc=rho1;
%     rho2_inc=rho2;
%
%     rho1_inc(i)=rho1_inc(i)+delta_rho;
%     rho2_inc(i)=rho2_inc(i)+delta_rho;
%
%     U1_c_inc=GravityPotentialTriDen93(x1',y1',z1',FV,xr,yr,zr,rho1_inc);
%     U2_c_inc=GravityPotentialTriDen93(x2',y2',z2',FV,xr,yr,zr,rho2_inc);
%
%     J1n(:,i)=(U1_c_inc-U1_c)'/delta_rho;
%     J2n(:,i)=(U2_c_inc-U2_c)'/delta_rho;
%
%     progressbar(i/Nelems);
%
% end
%
% progressbar(1);

% load jac3-1.mat;

%Jacobian numerical

tic;
J1=GravityPotentialTriDenDerFold93(x1',y1',z1',FV,xr,yr,zr,k);
J2=GravityPotentialTriDenDerFold93(x2',y2',z2',FV,xr,yr,zr,k);
t_deriv=toc


% tic;
% J1f=GravityPotentialTriDenDerFold93(x1',y1',z1',FV,xr,yr,zr,k);
% J2f=GravityPotentialTriDenDerFold93(x2',y2',z2',FV,xr,yr,zr,k);
% t_deriv=toc
% jacf=J1f-J2f;

% J2m=GravityPotentialTriDenDer(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);

jac=J1-J2;

clear J1 J2

% m=abs(J-jac);
%
% figure
% pcolor(log10(m));
% shading flat; colorbar;

% pcolor(log10(abs(J-jac)));shading flat
% sum(log10(abs(J(:)-jac(:)))<-8)

% Jb=GravityPotentialTriDenDer(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);

%% Inversion

L=U_o-U_c; % initial residuals

sumsq=sum(L.^2);

% cond=z1m>0;
% jac_t(:,1)=sum(jac(:,cond),2);
% jac_t(:,2)=sum(jac(:,~cond),2);
%
% jac_orig=jac;
% jac=jac_orig;

% reducing jacobian

% tic

kp=4^k;
% if (k>0)
%     jac_red=zeros(Nobs,Nelems/kp);
%     for n=1:Nelems/kp
%         jac_red(:,n)=sum(jac(:,n:Nelems/kp:Nelems),2);
%     end
%     jac=jac_red;
% end
% clear jac_red
%
% t_red_jac=toc

A=(jac'*jac);
b=jac'*L;

% drho1=A\b;

% drho1=((jac'*jac)^-1)*jac'*L;

range_crust=1900;

% Constrained inversion

G=sparse(3*Nelems/kp,Nelems/kp);

progressbar(0);

for i=1:Nelems/kp
    
    N = neighbors(TRi,i);
    
    G(3*i-2,i)=1;
    G(3*i-1,i)=1;
    G(3*i,  i)=1;
    
    G(3*i-2,N(1))=-1;
    G(3*i-1,N(2))=-1;
    G(3*i,  N(3))=-1;
    
    progressbar(i/(Nelems/kp));
end

progressbar(1)

lambda=0.01;

Astar=(A+lambda*(G'*G));

drho1=Astar\b;

std(drho1)

% lb=zeros(size(rho1))-range_crust;
% ub=zeros(size(rho1))+range_crust;
% drho1 = lsqlin(A,b,[],[],[],[],lb,ub);

% opts=optimset('TolFun',1e-16);
% [x,resnorm,residual,exitflag,output,lambda]  = lsqlin(A,b,[],[],[],[],lb,ub);
% std(x)


% rho1_s=rho1+drho1';

for n=1:Nelems/kp
    rho1_s(n:Nelems/kp:Nelems)=rho1(n:Nelems/kp:Nelems)+drho1(n);
end

rho2_s=ones(1,s(1))*rho_mantle-rho1_s;
rho3_s=ones(1,s(1))*rho_core-rho2_s-rho1_s;

% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho,refell);
U1_c=GravityPotentialTriDen93(x1',y1',z1',FV,xr,yr,zr,rho1_s);
U2_c=GravityPotentialTriDen93(x2',y2',z2',FV,xr,yr,zr,rho2_s);
U3_c=GravityPotentialTriDen93(x3',y3',z3',FV,xr,yr,zr,rho3_s);

U_c=U1_c+U2_c+U3_c;

clear U1_c U2_c U3_c

figure; hold on;
set(gca,'FontSize',20);
% plot((a_up-g_up)./g_up*100)
plot((L)./U_o*100,'.b');
plot((U_c-U_o)./U_o*100,'.r');
legend({'Initial residuals','Final residuals'},'FontSize',20);
ylabel('Potential residual [%] ','FontSize',20);

figure; hold on;
set(gca,'FontSize',20);
plot(rho1_s,'.');
plot(ones(size(rho1_s))*rho_crust+range_crust,'--k');
plot(ones(size(rho1_s))*rho_crust-range_crust,'--k');
plot(ones(size(rho1_s))*rho_crust,'-k');
ylabel('Density [kg/m^3]','FontSize',20);

% for fake generated potential
% figure; hold on; set(gca,'FontSize',20);
% plot(rho1_i,'bo');
% plot(rho1_s,'.r');
% ylabel('Density [kg/m^3]','FontSize',20);
% legend({'Initial density','Esimated density'},'FontSize',20);

% end

%% Plot Resuduals

[lambdare,fire]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );
[xre,yre,zre]=TriEllRadVec(fire,lambdare,aref,bref,cref,'xyz');

% Plotting initial residuals

% a1_up_re=RadGravityAccelerationTriDen(x1',y1',z1',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho1,refell);
% a2_up_re=RadGravityAccelerationTriDen(x2',y2',z2',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho2,refell);
% a3_up_re=RadGravityAccelerationTriDen(x3',y3',z3',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho3,refell);
%
% a_up_re=a1_up_re+a2_up_re+a3_up_re;
%
% [gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_obs,xre,yre,zre);
% [g_up_re,~,~]=GravityComponents(gx,gy,gz,xre,yre,zre,aref,cref);
%
% AGUaxes
% % pcolorm(fire,lambdare,g_up_re); shading interp;
% pcolorm(fire,lambdare,(g_up_re-a_up_re)*1e5); shading interp;
% % contourm(fir,lambdar,a,100,'Color','k'); shading interp;
% plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
% cbar=colorbar('FontSize',20);
% title('Residual anomaly','FontSize',20);
% ylabel(cbar,'Acceletation [mGal]','FontSize',20);
% % caxis(1e5*[min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])]);

% Plotting final residuals

% a_up_re=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho,refell);
U1_c_re=GravityPotentialTriDen93(x1',y1',z1',FV,xre,yre,zre,rho1_s);
U2_c_re=GravityPotentialTriDen93(x2',y2',z2',FV,xre,yre,zre,rho2_s);
U3_c_re=GravityPotentialTriDen93(x3',y3',z3',FV,xre,yre,zre,rho3_s);

U_c_re=U1_c_re+U2_c_re+U3_c_re;

clear U1_c_re U2_c_re U3_c_re

% [gx_re,gy_re,gz_re]=GravityAcceleration(mu,Rref,lmcosi_obs,xre,yre,zre);
% [g_up_re,g_east_re,g_north_re]=GravityComponents(gx_re,gy_re,gz_re,xre,yre,zre,aref,cref);

% U_o_re=GravityPotential2(mu,Rref,lmcosi_obs,xre,yre,zre);

[U_o_re,lambda_ell]=GravityPotentialEll4(mu,aref,RefEll,lmcosi,xre,yre,zre);

U_o_re=U_o_re/1000;

% [U_o_re,~,~,~]=plm2xyz(lmcosi_pot,fire(:)*180/pi,lambdare(:)*180/pi);
% U_o_re=reshape(U_o_re,size(fire));

%% Comparing initial potential (for fake generated initial potential)
% U1_o_re=GravityPotentialTriDen(x1',y1',z1',FV,xre,yre,zre,rho1_i);
% U2_o_re=GravityPotentialTriDen(x2',y2',z2',FV,xre,yre,zre,rho2_i);
% U3_o_re=GravityPotentialTriDen(x3',y3',z3',FV,xre,yre,zre,rho3_i);
%
% U_o_re=U1_o_re+U2_o_re+U3_o_re;

U1_init=GravityPotentialTriDen93(x1',y1',z1',FV,xre,yre,zre,rho1);
U2_init=GravityPotentialTriDen93(x2',y2',z2',FV,xre,yre,zre,rho2);
U3_init=GravityPotentialTriDen93(x3',y3',z3',FV,xre,yre,zre,rho3);

U_c_init=U1_init+U2_init+U3_init;

clear U1_init U2_init U3_init

U_res_init=U_o_re-U_c_init;

%% Initial Potential residuals
% [U_o_re2,lon,lat,Plm]=plm2xyz(lmcosi_pot,fire*180/pi,lambdare*180/pi);
AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,U_res_init); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
% plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Initial residuals','FontSize',20);
ylabel(cbar,'Potential [m^{2} s^{-2} ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
% caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

%% Plotting observed potential

AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,U_o_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from SH','FontSize',20);
ylabel(cbar,'Potential [m^2 s^-2 ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

%% Plotting calculated potential
AGUaxes
% pcolorm(fire,lambdare,a_up_re); shading interp;
pcolorm(fire,lambdare,U_c_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from polyhedra','FontSize',20);
ylabel(cbar,'Potential [m^2 s^-2 ]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

% res=(a_up_re-g_up_re)*1e5;
U_res_final=(U_c_re-U_o_re);

figure;
set(gca,'FontSize',20);
hist(U_res_final(:)./U_o_re(:)*100,20);
xlabel('Residual potential [%]','FontSize',20);

%% Final potential resuduals
AGUaxes
pcolorm(fire,lambdare,U_res_final); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
% plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Residual potential [m^{2} s^{-2} ]','FontSize',20);

%% Plot density map

PlotTriangles(xi,yi,zi,TRi.Triangulation,rho1_s)

%% SH expansion of density

MaxDegree=15;
Resolution=1;

lmcosi_den=xyz2plm(rho1_s,MaxDegree,'irr',fi1m'*180/pi,lambda1m'*180/pi);
[rho1_s_sh,lon_sh,lat_sh]=plm2xyz(lmcosi_den,Resolution);

[rho1_s_shm,lon_shm,lat_shm]=plm2xyz(lmcosi_den,fi1m*180/pi,lambda1m*180/pi);

[lon_sh,lat_sh]=meshgrid(lon_sh/180*pi,lat_sh/180*pi);

AGUaxes
pcolorm(lat_sh,lon_sh,rho1_s_sh); shading interp;
contourm(lat_sh,lon_sh,rho1_s_sh,'LevelStep',100,'ShowText','on','Color','k');

cbar=colorbar('FontSize',20);
ylabel(cbar,'Density [kg/m^{3}] ','FontSize',20);


%% Computing mass anomaly

V1=TetrahedraVolume(x1,y1,z1,FV);
V2=TetrahedraVolume(x2,y2,z2,FV);

Vcr=V1-V2;

Mcr=Vcr'.*rho1_s;

rho_crust_mean=sum(Mcr)/sum(Vcr);

Mcr=0;

Mcr=zeros(1,Nelems/kp);

for n=1:Nelems/kp
    Mcr(n)=sum(Vcr(n:Nelems/kp:Nelems)'.*rho1_s(n:Nelems/kp:Nelems));
end


PlotTriangles(x1o,y1o,z1o,FVo,Mcr)
PlotTriangles(x1o,y1o,z1o,FVo,rho1_s(1:Nelems/kp))

h_crust=mean(r1-r2)/1000;

%% Plot Vesta

figure; hold on; axis equal
set(gca,'FontSize',20);

% tri_fig=trisurf(FV,xi,yi,zi,ccol); shading interp
% set(tri_fig,'edgecolor','k','facecolor',ccol(1,:));
% trimesh_fig=trimesh(FV,xi,yi,zi,ri,'edgecolor','k','facecolor','none');

p = patch('Faces',FV,'Vertices',[x1/1000 y1/1000 z1/1000],'FaceColor','b',...
    'EdgeColor','none');

lighting phong;
light('Position',[-1 1 0],'Style','infinite');

clear cdata
cdata = rho1_s';

% cdata=(Vcr<0)*1
% n=317;
% cdata = (jac(:,n)<0)*1;

% cdata = (rho1_inc+rho2_inc)';
set(p,'FaceColor','flat',...
    'FaceVertexCData',cdata,...
    'CDataMapping','scaled')

% caxis([2000 4000]);

axis equal

cbar=colorbar('FontSize',20);
ylabel(cbar,'Density [kg/m^3] ','FontSize',20);

xlabel('x [km]','FontSize',20);
ylabel('y [km]','FontSize',20);
zlabel('z [km]','FontSize',20);

box on;

%% ShowVesta

% [rho1_s_sh_hr,~,~]=plm2xyz(lmcosi_den,0.05);
% rho1_s_sh_hr=circshift(rho1_s_sh_hr,[0 3600]);
% set(h,'CData',rho1_s_sh_hr');

