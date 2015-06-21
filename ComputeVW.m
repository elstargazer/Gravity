function [V,W]=ComputeVW(x,y,z,Rref,n,m)

r=sqrt(x.*x+y.*y+z.*z);
factor=(Rref./r)^(n+1);

lambda=atan2(y,x);


Pnm=legendre(n,z./r);

Pnm=Pnm(m+1);

V=factor.*Pnm.*cos(m*lambda)*(-1)^m;
W=factor.*Pnm.*sin(m*lambda)*(-1)^m;