function lmcosi_sub = FindSubRelief(lmcosi_g,lmcosi_t,GM,Rref,rho1,rho2,r2,T)
% This is a programs that computes subsurface interface given
% observed gravity/shape & assumed densities.

%% Parameters 
G = 6.67384e-11;
M = GM/G;

hmaxt  = 4;      % max power of topography
nmaxg  = 2;      % max degree of gravity
nmaxgt = nmaxg;  % max degree of gravity from shape
nmaxt  = 10;     % max degree of shape
res    = 1;      % resolution [deg]
nc     = 2;      % critical degree
rtol   = 1;      % tolerance for subsurface relief [m]
aref   = 481000; % reference for Bouguer anomaly
cref   = 446000; % reference for Bouguer anomaly

%% Computing hydrostatic state

[ri,lon,lat] = plm2xyz(lmcosi_t,res);
[lon,lat]    = meshgrid(lon,lat);
lon = lon/180*pi;
lat = lat/180*pi;

r1 = lmcosi_t(1,3);

[fh,~]=HydrostaticStateExact2l(r1,r2,T,rho1,rho2,0.1,0.1);

fp1 = fh(1);
fp2 = fh(2);

drho      = rho2 - rho1;
[a2,~,c2] = fr2abc(r2,fp2,0);
M2        = muEllipsoid(a2,a2,c2,drho)/G;
M1        = M - M2;

%% Computing Bouguer anomaly

ri2 = TriEllRadVec(lat,lon,a2,a2,c2,'rad');

lmcosi_t2_ell  = xyz2plm(ri2,nmaxt);

lmcosi_gt2 = SHRotationalEllipsoid(a2,c2,nmaxg,Rref); 
lmcosi_gt1 = Topo2Grav(ri,Rref,nmaxt,nmaxgt,hmaxt);

lmcosi_gt        = lmcosi_gt1;
lmcosi_gt(:,3:4) = (M1*lmcosi_gt1(:,3:4)+M2*lmcosi_gt2(:,3:4))/M;

lmcosi_ba        = lmcosi_g;
lmcosi_ba(:,3:4) = lmcosi_ba(:,3:4) - lmcosi_gt(:,3:4);
lmcosi_ba(4,3) = 0;
%% Computing Bouguer anomaly

% [xref,yref,zref] = TriEllRadVec(lat,lon,aref,aref,cref,'xyz');
% 
% [ax_ba,ay_ba,az_ba]=GravityAcceleration(...
%     GM,Rref,lmcosi_ba,xref,yref,zref);
% 
% [g_up_ba,~,~]=GravityComponents(...
%     ax_ba,ay_ba,az_ba,xref,yref,zref,aref,cref);
%     
% AGUaxes;
% pcolorm(lat*180/pi,lon*180/pi,(g_up_ba)*1e5);
% cbar = colorbar('FontSize',20);
% ylabel(cbar,'Bouguer anomaly [mGal]','FontSize',20);


%% Computing subrelief

n = lmcosi_ba(:,1);
A = M*(2.*n+1).*((Rref./r2).^n)./...
    (4*pi*drho.*(r2^2));
A = [A A];
ll = 1./(M*(2.*nc+1).*((Rref./r2).^nc)./...
    (4*pi*drho.*(r2^2))).^2;
w = (1+ll.*(A.^2)).^(-1);

lmcosi_sub = CreateEmptylmcosi(nmaxg);
lmcosi_sub(:,3:4) = w.*lmcosi_ba(:,3:4).*A;
lmcosi_sub_corr = lmcosi_sub;

dro = rtol + 1;

ri2_old = ri2 - r2;

while (dro > rtol)
    
    lmcosi_sub_corr(:,3:4) = 0;
    
    for h=2:hmaxt
        
        ri2 = plm2xyz(lmcosi_sub,res);
       
        lmcosi_sub_h = xyz2plm(ri2.^h,nmaxgt);
        
        P = ones(numel(n),1);
        for j=1:h
            P=P.*(n+4-j);
        end
        
        P = [P P];        
        B = P./((r2.^h).*factorial(h).*([n n]+3));
            
        lmcosi_sub_corr(:,3:4) = lmcosi_sub_corr(:,3:4) + ...
            w.*r2.*lmcosi_sub_h(:,3:4).*B;
        lmcosi_sub_corr(1,3) = 0;
    end
    
    lmcosi_sub(:,3:4) = lmcosi_sub(:,3:4) - lmcosi_sub_corr(:,3:4);
%     lmcosi_sub(1,3) = r2;
    
    
    ri2 = plm2xyz(lmcosi_sub,res);
    
    dro = ri2 - ri2_old;
    dro = max(abs(dro(:)));
    dro
    
    ri2_old = ri2;
end

lmcosi_sub(:,3:4) = lmcosi_t2_ell(:,3:4) + lmcosi_sub(:,3:4);



















