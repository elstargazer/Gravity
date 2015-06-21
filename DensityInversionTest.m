ccc

%% Load gravity model

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

N_trunc=15;

lmcosi_obs=TruncateGravityModel(lmcosi_obs_full,N_trunc,0);



%% Load topography model

NTess=3;
filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
[FV,xi,yi,zi]=Grid2TriShapeModel(filename,NTess);

s=size(FV);

xi=xi*1000;
yi=yi*1000;
zi=zi*1000;

xim=(xi(FV(:,1))+xi(FV(:,2))+xi(FV(:,3)))/4;
yim=(yi(FV(:,1))+yi(FV(:,2))+yi(FV(:,3)))/4;
zim=(zi(FV(:,1))+zi(FV(:,2))+zi(FV(:,3)))/4;

[lambdai,fii,ri]=cart2sph(xi,yi,zi);
[lambdaim,fiim,rim]=cart2sph(xim,yim,zim);

rho_mean=3457;
rho=ones(1,s(1))*rho_mean;

%% Create reference surface

aref=300000;
cref=300000;
TR=IcosahedronMesh;

TR_2r=SubdivideSphericalMesh(TR,NTess+2); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FVr=TR_2r.Triangulation;

xr=aref*TR_2r.X(:,1);
yr=aref*TR_2r.X(:,2);
zr=cref*TR_2r.X(:,3);

xr=(xr(FVr(:,1))+xr(FVr(:,2))+xr(FVr(:,3)))/3;
yr=(yr(FVr(:,1))+yr(FVr(:,2))+yr(FVr(:,3)))/3;
zr=(zr(FVr(:,1))+zr(FVr(:,2))+zr(FVr(:,3)))/3;


[lambdar,fir,~]=cart2sph(xr,yr,zr);

[xr,yr,zr]=sph2cart(lambdar,fir,1);

xr=aref*xr;
yr=aref*yr;
zr=cref*zr;

[~,~,rr]=cart2sph(xr,yr,zr);

lmcosi_pot=plm2pot(AddZeroHarm(lmcosi_obs,1),aref,mu,Rref,1,'nothing');

%% Plot

figure; hold on; axis equal

tri_fig=trisurf(FV,xi,yi,zi,ri); shading interp
set(tri_fig,'edgecolor','k'); 
trimesh_fig=trimesh(FV,xi,yi,zi,ri,'edgecolor','k','facecolor','none');

% surf(xr,yr,zr,a); 

plot3(xr,yr,zr,'.k','MarkerSize',10);
axis(1.4*[-300000 300000 -300000 300000 -300000 300000])
alpha(0.2);

%% Computing gravity
% [gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_obs,xr,yr,zr);
% [g_up,g_east,g_north]=GravityComponents(gx,gy,gz,xr,yr,zr,aref,cref);

% U_o=GravityPotential(mu,Rref,lmcosi_obs,xr,yr,zr);

[U_o,~,~,~]=plm2xyz(lmcosi_pot,fir*180/pi,lambdar*180/pi);



% gtot=sqrt(gx.*gx+gy.*gy+gz.*gz);

% [ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% atot=sqrt(ax.*ax+ay.*ay+az.*az);

% U_c=GravityPotentialTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);


% d2=DeltaSqGravity(rho,xi,yi,zi,FV,xr,yr,zr,gx,gy,gz)

%% Numerical Jacobian

% Jacobian linear
% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
% a_up0=a_up;

U_c=GravityPotentialTri93(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_mean);
U_c0=U_c;


% jac=zeros(numel(xr),numel(rho));

Nelems=numel(rho);


% rho_inc=rho;
% delta_rho=1000;
% progressbar(0);
% 
% 
% 
% for i=1:Nelems   
%     
%     rho_inc=rho;
%     rho_inc(i)=rho_inc(i)+delta_rho;    
% %     a_up1=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_inc); 
%     U_c1=GravityPotentialTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_inc); 
% %     jac(:,i)=(a_up1-a_up0)'/delta_rho;   
%     jac(:,i)=(U_c1-U_c)'/delta_rho;   
%     progressbar(i/Nelems);  
%     
% end
% 
% progressbar(1);


jac=GravityPotentialTriDenDer93(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);


% J=GravityPotentialTriDenDer(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);

% m=abs(J-jac);
% 
% figure
% pcolor(log10(m));
% shading flat; colorbar;
% 
% find(log10(m(1,:))<-8)

% Jb=GravityPotentialTriDenDer(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr);

% jacobian numerical

% fun = @(rho_s) RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_s);
% 
% % d=distance2(xr,yr,zr,xim,yim,zim);
% % for i=1:5
% 
% 
% tic
% [jac,err] = jacobianest(fun,rho);
% toc

% jac(d>150000)=0;

% jac_t=jac;
% jac_t(abs(jac)>3e-7)=0;

% L=g_up-a_up;
L=U_o-U_c;

sumsq=sum(L.^2);

A=(jac'*jac);
b=jac'*L;

drho=A\b;

% lb=zeros(size(rho))-1000;
% ub=zeros(size(rho))+1000;% 
% drho = lsqlin(A,b,[],[],[],[],lb,ub);

rho_s=rho+drho';


% a_up=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);
U_c=GravityPotentialTriDen93(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho_s);


figure; hold on;
% plot((a_up-g_up)./g_up*100)
plot((L)./U_o*100,'.b')
plot((U_c-U_o)./U_o*100,'.r')

figure
plot(rho_s)
ylabel('Density [kg/m^3]','FontSize',20);

% end

%% Plot Resuduals

FiStep=2;
LambdaStep=2;
[lambdare,fire]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );

% aref=293200;
% cref=266500;

xre=aref*cos(lambdare).*cos(fire);
yre=aref*sin(lambdare).*cos(fire);
zre=cref*sin(fire);

rre=sqrt(xre.*xre+yre.*yre+zre.*zre);

% a_up_re=RadGravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho);
U_c_re=GravityPotentialTriDen93(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho_s);
U_c_init=GravityPotentialTri93(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho_mean);
% [gx_re,gy_re,gz_re]=GravityAcceleration(mu,Rref,lmcosi_obs,xre,yre,zre);
% [g_up_re,g_east_re,g_north_re]=GravityComponents(gx_re,gy_re,gz_re,xre,yre,zre,aref,cref);
% U_o_re=GravityPotential(mu,Rref,lmcosi_obs,xre,yre,zre);

[U_o_re,~,~,~]=plm2xyz(lmcosi_pot,fire(:)*180/pi,lambdare(:)*180/pi);
U_o_re=reshape(U_o_re,size(fire));


AGUaxes
% pcolorm(fire,lambdare,g_up_re); shading interp;
pcolorm(fire,lambdare,U_o_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from SH','FontSize',20);
ylabel(cbar,'Acceletation [m/s^2]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])

AGUaxes
% pcolorm(fire,lambdare,a_up_re); shading interp;
pcolorm(fire,lambdare,U_c_re); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
title('Gravity from polyhedra','FontSize',20);
ylabel(cbar,'Acceletation [m/s^2]','FontSize',20);
% caxis([min([g_up_re(:); a_up_re(:)]) max([g_up_re(:); a_up_re(:)])])
caxis([min([U_o_re(:); U_c_re(:)]) max([U_o_re(:); U_c_re(:)])])


% res=(a_up_re-g_up_re)*1e5;
res_init=(U_c_init-U_o_re);
res=(U_c_re-U_o_re);

figure; 
set(gca,'FontSize',20);
hist(res(:)./U_o_re(:)*100,20);

AGUaxes
pcolorm(fire,lambdare,res_init); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Initial residual potential [mGal]','FontSize',20);


AGUaxes
pcolorm(fire,lambdare,res); shading interp;
% contourm(fir,lambdar,a,100,'Color','k'); shading interp;
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',9);
cbar=colorbar('FontSize',20);
ylabel(cbar,'Final residual potential [mGal]','FontSize',20);

%% Plot Vesta

figure; hold on; axis equal

% tri_fig=trisurf(FV,xi,yi,zi,ccol); shading interp
% set(tri_fig,'edgecolor','k','facecolor',ccol(1,:)); 
% trimesh_fig=trimesh(FV,xi,yi,zi,ri,'edgecolor','k','facecolor','none');

p = patch('Faces',FV,'Vertices',[xi yi zi],'FaceColor','b');

clear cdata
cdata = rho_s';
set(p,'FaceColor','flat',...
'FaceVertexCData',cdata,...
'CDataMapping','scaled')

caxis([2000 4000]);

axis equal

colorbar

% rho(abs((rho-rho_mean))>500)=rho_mean










