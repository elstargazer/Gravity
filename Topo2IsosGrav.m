function lmcosi_isosgrav = Topo2IsosGrav(...
    ri,Rref,D,rho_crust,rho_mantle,rho_mean,nmaxt,nmaxgt,hmaxt)


% conversion from topograhy coefficients to isotatically compensated
% gravity coefficients

delta_rho = rho_mantle - rho_crust;

tic
progressbar('Conversion to isostatic gravitational harmonics');

% wrap up the shape
ri(:,1)=ri(:,end);

% expand shape in spherical harmonics
[lmcosi_topo,~]=xyz2plm(ri,nmaxt,'im',[],[],[]);
n=lmcosi_topo(:,1);
R0=lmcosi_topo(1,3);

Coef1=3./(2*n+1);

lmcosi_grav=CreateEmptylmcosi(nmaxt);
lmcosi_isos=CreateEmptylmcosi(nmaxt);

% normalize shape coefficients
lmcosi_grav(:,3)=lmcosi_topo(:,3)/R0;
lmcosi_grav(:,4)=lmcosi_topo(:,4)/R0;

IsosCoef=((R0-D)/R0).^n;

lmcosi_isos(:,3)=lmcosi_topo(:,3)/R0.*IsosCoef;
lmcosi_isos(:,4)=lmcosi_topo(:,4)/R0.*IsosCoef;

% subtract mean radius
ri=ri-R0;

progressbar(1/hmaxt);

for h=2:hmaxt
    
    P=ones(numel(n),1);
    
    for j=1:h
        P=P.*(n+4-j);
    end
    
    Coef2=P./((R0^h)*factorial(h).*(n+3));
    rin=ri.^h;
    [lmcosi_topo_power,~]=xyz2plm(rin,nmaxt,'im',[],[],[]);
    
    lmcosi_grav(:,3)=lmcosi_grav(:,3)+Coef2.*lmcosi_topo_power(:,3);
    lmcosi_grav(:,4)=lmcosi_grav(:,4)+Coef2.*lmcosi_topo_power(:,4);
    
    IsosCoef = (-1)^(h-1)*(((R0-D)/R0).^(n - (3*h-1))).*(rho_crust/delta_rho).^(h-1);
    
    lmcosi_isos(:,3)=lmcosi_isos(:,3)+Coef2.*IsosCoef.*lmcosi_topo_power(:,3);
    lmcosi_isos(:,4)=lmcosi_isos(:,4)+Coef2.*IsosCoef.*lmcosi_topo_power(:,4);
    
    progressbar(h/hmaxt);
end

Coef3=(R0/Rref).^n;

lmcosi_grav(:,3)=lmcosi_grav(:,3).*Coef1.*Coef3;
lmcosi_grav(:,4)=lmcosi_grav(:,4).*Coef1.*Coef3;

lmcosi_isos(:,3)=lmcosi_isos(:,3).*Coef1.*Coef3;
lmcosi_isos(:,4)=lmcosi_isos(:,4).*Coef1.*Coef3;

lmcosi_isosgrav = lmcosi_grav;

lmcosi_isosgrav(:,3) = lmcosi_isosgrav(:,3) - lmcosi_isos(:,3);
lmcosi_isosgrav(:,4) = lmcosi_isosgrav(:,4) - lmcosi_isos(:,4);

lmcosi_isosgrav(1,3)=1;
lmcosi_isosgrav=TruncateGravityModel(lmcosi_isosgrav,nmaxgt,1);

toc
