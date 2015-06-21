ccc
%% Load gravity model
EllFileName='/Users/antonermakov/Dawn/Gravity/VESTAELL/nomfil_VESEH20D3.txt';

G=6.67384e-11;

Ntrunc=17;
NTess=5;

lmcosi=ReadEllHarmModel(EllFileName);

% Plot spectum
[l,sdl,sdl_error]=EllipsoidHarmonicSpectrum(lmcosi);

figure; hold on;
set(gca,'FontSize',20);

plot(l,sdl,'-o','LineWidth',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(l,sdl_error,'r-o','LineWidth',3,'MarkerFaceColor','k','MarkerEdgeColor','k');

set(gca,'YScale','log');

xlabel('Degree','FontSize',20);
ylabel('Spectral Density','FontSize',20);
box on;

lmcosi=TruncateEllModel(lmcosi,Ntrunc);

% nc=8;
% pc=5;
% % selind=(lmcosi(:,1)==nc) 
% selind=(lmcosi(:,1)==nc) & (lmcosi(:,2)==pc);
% lmcosi(~selind,3)=0;
% lmcosi(selind,3)=1;

a=303.860;
b=288.660;
c=246.460;

% Re=265.000;
Re=a;

% a=3;
% b=2;
% c=1;

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

RefEll=[a b c];

mu=17.28830802930000;

angle_fi=1.209/180*pi;
angle_theta=-0.558/180*pi;
angle_psi=129.209/180*pi;

S=rot(angle_psi,3)*rot(angle_theta,2)*rot(angle_fi,1);

%% Generate surface

% x=400000;
% y=400000;
% z=400000;

% aref=a;
% bref=b;
% cref=c;

% aref=1;
% bref=1;
% cref=1;

Shift=0.1;
GridStep=2;

fi=[(-90+Shift:GridStep:90) 90-Shift];
lambda=[Shift:GridStep:360 360-Shift];

% fi=-90:GridStep:90;
% lambda=0:GridStep:360;

[fii,lambdai]=meshgrid(fi/180*pi,lambda/180*pi);

Npts=numel(fii);

[x,y,z]=TriEllRadVec(fii,lambdai,a,b,c);
r=sqrt(x.*x+y.*y+z.*z);

U0=mu./r;
U0mean=mean(U0(:))*1e6;


%% Plotting ellipsoid
AGUaxes
pcolorm(fii,lambdai,r)
shading interp
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius-Vector Magnitude [km]','FontSize',20);


[xrotell,yrotell,zrotell]=RotateByMatrix(x,y,z,S');
[lambdairotell,fiirotell,rrotell]=cart2sph(xrotell,yrotell,zrotell);

AGUaxes
set(gca,'FontSize',20);
pcolorm(fiirotell,lambdairotell,rrotell)
shading interp
cbar=colorbar('FontSize',20);
ylabel(cbar,'Radius-Vector Magnitude [km] ','FontSize',20);

%% Computing potential

tic
[U1,lambda_ell]=GravityPotentialEll6(mu,Re,RefEll,lmcosi,x,y,z); U1=U1*1e6;
toc

tic
[U,lambda_ell]=GravityPotentialEll4(mu,Re,RefEll,lmcosi,x,y,z); U=U*1e6;
toc


AGUaxes
set(gca,'FontSize',20);
pcolorm(fiirotell,lambdairotell,U1)
shading interp
cbar=colorbar('FontSize',20);
ylabel(cbar,'Potential [m^2 s^{-2}] ','FontSize',20);



%% Reading Covariance

CovFileName='/Users/antonermakov/Dawn/Gravity/VESTAELL/cov_VESEH20D3.txt';
Cov=ReadEllCovariance(CovFileName);

Cov(EllIndex(Ntrunc+1,1):end,:)=[];
Cov(:,EllIndex(Ntrunc+1,1):end)=[];

%% Propagating covariance
tic
[sigma,U2_t,var_t]=GravityPotentialEllCov(mu,Re,RefEll,lmcosi,Cov,x,y,z); 
% sigma=sigma*1e6;
% var_t=var_t*1e6;
% U2_t=U2_t*1e6;

toc

AGUaxes
set(gca,'FontSize',20);
pcolorm(fiirotell,lambdairotell,sigma)
shading interp
cbar=colorbar('FontSize',20);
ylabel(cbar,'Potential [m^2 s^{-2}] ','FontSize',20);

% tic
% [ax,ay,az,lambda_ell]=GravityAccelerationEll(mu,Re,RefEll,lmcosi,x,y,z);
% [au,ae,an]=GravityComponents(-1000*ax,-1000*ay,-1000*az,x,y,z,1,1);
% at=1000*sqrt(ax.*ax+ay.*ay+az.*az);
% toc;


%% Computing degree strength

deg_str=zeros(Npts,1);

U2_t=U2_t*1e6;

sig_pow=(U2_t(:,3:end));
sig_pow_mean=(mean(sig_pow));


deg=2:Ntrunc;
logdeg=log10(deg);

% figure; hold on;

progressbar('Computing degree strength');

for i=1:Npts

    p_s=polyfit(logdeg,log10(U2_t(i,3:end)),1);    
    U_t_fit=polyval(p_s,logdeg);
    
%       plot(logdeg,log10(sig_pow(i,:)),'b','LineWidth',3);   
%       plot(logdeg,log10(var_t(i,3:end)),'r','LineWidth',3);    
%       plot(logdeg,U_t_fit,'b--','LineWidth',3);
%      ylim([-5 5]);
    
    [deg_str(i),~,~,~] = intersections(logdeg,log10(var_t(i,3:end)),...
                                      logdeg,U_t_fit);   
%     [deg_str(i),~,~,~] = intersections(logdeg,log10(sigma_t(i,3:end).^2),...
%                                        logdeg,log10(sig_pow_mean));                              
%     plot(n_intesect,y0,'*r','MarkerSize',30); 
     progressbar(i/Npts);
% 
%     drawnow;
%     unplot(3)
end

progressbar(1);

deg_str=reshape(deg_str,size(fii));
deg_str=10.^deg_str;

AGUaxes
pcolorm(fii,lambdai,sigma)
shading interp;
cbar=colorbar('FontSize',20);

AGUaxes
pcolorm(fii,lambdai,U)
shading interp;
cbar=colorbar('FontSize',20);

AGUaxes
pcolorm(fii,lambdai,deg_str)
shading interp;
cbar=colorbar('FontSize',20);

AGUaxes
pcolorm(fii,lambdai,r)
shading interp;
cbar=colorbar('FontSize',20);

%% Movie

% vidObj = VideoWriter('EllModelError.avi');
% open(vidObj);
% 
% AGUaxes
% plot1=pcolorm(fii,lambdai,reshape(sigma_t(:,1),size(fii)));
% shading interp;
% cbar=colorbar('FontSize',20);
% title('Degree = 0 ','FontSize',20);
% 
% ylabel(cbar,'Potential Error [m^{2} s^{-2}] ','FontSize',20);
% 
% currFrame = getframe(gcf);
% writeVideo(vidObj,currFrame);
% 
% for n=2:Ntrunc+1    
%     set(plot1,'CData',reshape(sigma_t(:,n),size(fii)));
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
%     title(['Degree = ' num2str(n-1) ' '],'FontSize',20);
% end
% 
% close(vidObj);


% U=U*1e6;

% n=2;
% p=5;
% Ya=SurfaceEllipsoidalHarmonic2(p,RefEll,x,y,z);
% Yn=SurfaceEllipsoidalHarmonic(n,p,RefEll,x,y,z);
% Umean=mean(U(:));
% Umean/U0mean
% 
% eps=1e-5;
% n=7;
% p=1;
% 
% [P,Pderiv] = calcPM(lambda_ell(:,2), n, a, b, c);
% [P2,~] = calcPN(lambda_ell(:,2)+eps, n, a, b, c);
% 
% [d11,d12,d13,d21,d22,d23,d31,d32,d33]=EllRectDer(lambda_ell,x,y,z,a,b,c);
% 
% % 
% % Pderiv2=(P2-P)/eps;
% 
% 
% lambda_ell = approxCartToEllnosign(a, b, c, [x(:) y(:) z(:)]);   
% lambda_ell_xplus = approxCartToEllnosign(a, b, c, [x(:)+eps y(:) z(:)]);   
% lambda_ell_yplus = approxCartToEllnosign(a, b, c, [x(:) y(:)+eps z(:)]);   
% lambda_ell_zplus = approxCartToEllnosign(a, b, c, [x(:) y(:) z(:)+eps]);   
% 
% dlambda_ellx=(lambda_ell_xplus-lambda_ell)./eps;
% dlambda_elly=(lambda_ell_yplus-lambda_ell)./eps;
% dlambda_ellz=(lambda_ell_zplus-lambda_ell)./eps;

%% Plotting potential
% computed potential

% pcolorm(fii,lambdai,reshape((Pderiv2(:,p)-Pderiv(:,p))./Pderiv(:,p),size(x)))
% pcolorm(fii,lambdai,reshape(P(:,1),size(x)));
% pcolorm(fii,lambdai,reshape(d11,size(x)));

% d11d=(dlambda_ellx(:,1)-d11)./dlambda_ellx(:,1);
% d12d=(dlambda_ellx(:,2)-d12)./dlambda_ellx(:,2);
% d13d=(dlambda_ellx(:,3)-d13)./dlambda_ellx(:,3);
% 
% d21d=(dlambda_elly(:,1)-d21)./dlambda_elly(:,1);
% d22d=(dlambda_elly(:,2)-d22)./dlambda_elly(:,2);
% d23d=(dlambda_elly(:,3)-d23)./dlambda_elly(:,3);
% 
% d31d=(dlambda_ellz(:,1)-d31)./dlambda_ellz(:,1);
% d32d=(dlambda_ellz(:,2)-d32)./dlambda_ellz(:,2);
% d33d=(dlambda_ellz(:,3)-d33)./dlambda_ellz(:,3);
% 
% 
% AGUaxes
% pcolorm(fii,lambdai,reshape(d33d,size(x)));
% shading interp;
% cbar=colorbar('FontSize',20);

% caxis([-1e-3 1e-3])


% d33d(abs(d33d)==Inf)=[];
% hist(log10(abs(d33d)),500);



%% Shape model
%mantle
% a2=257.000;
% c2=207.000;
% 
% % core
% a3=117.000;
% c3=105.000;
% 
% rhom=3457;
% rho_core=7400;
% rho_mantle=3160;
% rho_crust=2885*1.037231332824880;
% 
% filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
% 
% [TR,x1,y1,z1]=Grid2TriShapeModel(filename,NTess);
% FV=TR.Triangulation;
% s=size(FV);
% 
% 
% rho1=rho_crust;
% rho2=rho_mantle-rho1;
% rho3=rho_core-rho2-rho1;
% 
% [lambda1,fi1,r1]=cart2sph(x1,y1,z1);
% 
% [x2,y2,z2]=TwoEllRadVec(lambda1,fi1,a2,c2);
% [x3,y3,z3]=TwoEllRadVec(lambda1,fi1,a3,c3);
% 
% % dx=-3.333695212899673e+02/1000;
% % dy=-1.407808477824721e+03/1000;
% % dz=4.442896381860520/1000;
% 
% dx=0;
% dy=0;
% dz=0;
% 
% x1=x1-dx;
% y1=y1-dy;
% z1=z1-dz;
% 
% x2=x2-dx;
% y2=y2-dy;
% z2=z2-dz;
% 
% x3=x3-dx;
% y3=y3-dy;
% z3=z3-dz;
% 
% r1=sqrt(x1.*x1+y1.*y1+z1.*z1);
% r2=sqrt(x2.*x2+y2.*y2+z2.*z2);
% r3=sqrt(x3.*x3+y3.*y3+z3.*z3);
% 
% 
% % rp=TriEllRadVec(fii,lambdai,293,293,266);
% % [xp,yp,zp]=sph2cart(lambdai,fii,rp);
% 

% 
% [x1r,y1r,z1r]=RotateByMatrix(x1,y1,z1,S);
% [x2r,y2r,z2r]=RotateByMatrix(x2,y2,z2,S);
% [x3r,y3r,z3r]=RotateByMatrix(x3,y3,z3,S);
% 
% tic
% [U1_c,a1x_c,a1y_c,a1z_c]=GravityPotAccTri(1000*x1r',1000*y1r',1000*z1r',FV,1000*x,1000*y,1000*z,rho1);
% [U2_c,a2x_c,a2y_c,a2z_c]=GravityPotAccTri(1000*x2r',1000*y2r',1000*z2r',FV,1000*x,1000*y,1000*z,rho2);
% [U3_c,a3x_c,a3y_c,a3z_c]=GravityPotAccTri(1000*x3r',1000*y3r',1000*z3r',FV,1000*x,1000*y,1000*z,rho3);
% t_pot_calc=toc
% 
% U_c=(U1_c+U2_c+U3_c);
% 
% ax_c=a1x_c+a2x_c+a3x_c;
% ay_c=a1y_c+a2y_c+a3y_c;
% az_c=a1z_c+a2z_c+a3z_c;
% 
% [au_c,ae_c,an_c]=GravityComponents(ax_c,ay_c,az_c,x,y,z,1,1);
% 
% at_c=sqrt(ax_c.*ax_c+ay_c.*ay_c+az_c.*az_c);
% 
% 
% % U=U+2*(mean(Ut(:))-mean(U(:)));
% 
% AGUaxes
% pcolorm(fii,lambdai,(au_c)*1e5);
% shading interp;
% cbar=colorbar('FontSize',20);
% 
% AGUaxes
% pcolorm(fii,lambdai,(U-U_c));
% shading interp;
% cbar=colorbar('FontSize',20);
% 
% MaxDegree=80;
% lmcosi_shape=xyz2plm(r1,MaxDegree,'irr',fi1*180/pi,lambda1*180/pi);
% r1_sh=plm2xyz(lmcosi_shape,fii(:)*180/pi,lambdai(:)*180/pi);
% r1_sh=reshape(r1_sh,size(fii));
% 
% AGUaxes
% pcolorm(fii,lambdai,((r1_sh-r)<0)*1);
% shading interp;
% cbar=colorbar('FontSize',20);
 
% U_c_mean=mean(U_c(:))

%% Movie harmonics
% vidObj = VideoWriter('EllipsoidalHarmonics.avi');
% open(vidObj);
% 
% gamma = computeRomainNormalizationConstants(15, a, b, c);
% 
% AGUaxes
% 
% for i=1:size(lmcosi,1)
%     
%     n=lmcosi(i,1);
%     p=lmcosi(i,2);
%     class=LameClass(n,p);
% 
%     pcolorm(fii,lambdai,reshape(E_mu(:,i).*E_nu(:,i),size(U))./sqrt(gamma(EllIndex(n,p))));
%     shading interp;
%     cbar=colorbar('FontSize',20);
%    
%     title(['{n,p} = ' num2str([n p]) '; Class = ' class ' '],'FontSize',20);
%    
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
%     unplot
% end
% 
% 
% close(vidObj);

% AGUaxes
% pcolorm(fii,lambdai,Ya);
% shading interp;
% cbar=colorbar('FontSize',20);
% 
% 
% AGUaxes
% pcolorm(fii,lambdai,Ya./Yn);
% shading interp;
% cbar=colorbar('FontSize',20);

% ylabel(cbar,'Potential [m^{2} s^{-2}]','FontSize',20);


% point mass potential
% AGUaxes
% pcolorm(fii,lambdai,U0*1e6);
% shading interp;
% cbar=colorbar('FontSize',20);
% ylabel(cbar,'Potential [m^{2} s^{-2}]','FontSize',20);


%% Plot Coordinates
% 
% lambda1=lambda_ell(:,1);
% lambda2=lambda_ell(:,2);
% lambda3=lambda_ell(:,3);
% 
% AGUaxes
% pcolorm(fii,lambdai,reshape(E_nu(:,EllIndex(3,5)),size(U))); shading interp
% colorbar('FontSize',20);


% figure; hold on;
% set(gca,'FontSize',20);
% pcolor(reshape(lambda_ell(:,2),size(U)),...
%        reshape(lambda_ell(:,3),size(U)),...
%        reshape(E_nu(:,EllIndex(1,1)),size(U))); shading interp
% 
% xlabel('\mu','FontSize',20);
% ylabel('\nu','FontSize',20);

















 
