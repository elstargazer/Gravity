ccc;

filename_t = '~/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150702_GRAVITY_SPC/tsh_SPC150708a_n20.nml';
filename_gt = '~/Dawn/CeresGravityModel/gsh_SPC150708a_n20.nml';
[lmcosi_gtr,Rref] = ReadRyanSH(filename_gt);
[lmcosi_tr,Rref] = ReadRyanSH(filename_t);
lmcosi_tr(:,3:4) = lmcosi_tr(:,3:4) * Rref; 

GM    = 62.6253e9;
lmcosi_gtr(1,3) = 1;

figure; hold on; box on;
set(gca,'YScale','log');
set(gca,'XScale','log')
set(gca,'FontSize',20);
ylabel('Power','FontSize',20);
xlabel('Degree','FontSize',20);

[sdl1,l1]=plm2spec(lmcosi_gt1);
[sdlr,lr]=plm2spec(lmcosi_gtr);

plot(l1,sdl1,'r-');
plot(lr,sdlr,'b-');

[n,Z] = SphericalHarmonicAdmittance(lmcosi_gtr,lmcosi_tr,GM,Rref);
r = SphericalHarmonicCorrelation(lmcosi_tr,lmcosi_gtr);

aref  = 481000;
cref  = 446000;

lat = (-90:1:90)/180*pi;
lon = (-180:1:180)/180*pi;
[lat,lon]=meshgrid(lat,lon);   
[xref,yref,zref]=TriEllRadVec(lat,lon,aref,aref,cref,'xyz');

[ax,ay,az]=GravityAcceleration(...
    GM,Rref,lmcosi_gtr,xref,yref,zref);
[g_upr,g_east,g_north]=GravityComponents(...
    ax,ay,az,xref,yref,zref,aref,cref);

AGUaxes;
pcolorm(lat*180/pi,lon*180/pi,g_up);


[ax,ay,az]=GravityAcceleration(...
    GM,Rref,lmcosi_gt1,xref,yref,zref);
[g_up,g_east,g_north]=GravityComponents(...
    ax,ay,az,xref,yref,zref,aref,cref);

AGUaxes;
pcolorm(lat*180/pi,lon*180/pi,(g_up-g_upr)*1e5);
colorbar





