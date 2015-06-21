function [ell_estim,rms]=FindRotationalEllipsoidRaw(r,fi,lambda,ell)

options=optimset;
options=optimset(options,'UseParallel','always');

[ell_estim,rms]=fminsearch(@(ell) FindParameter2(ell, r,fi,lambda),ell,options);


function SqDist=FindParameter2(ell, r,fi,lambda)

[x,y,z]=sph2cart(lambda,fi,r);

e=sqrt(1-(ell(2)^2)/(ell(1)^2));


[B,L,H]=XYZ2BLH(x,y,z,ell(1),e);

SqDist=sum(sum(H.^2));
ell