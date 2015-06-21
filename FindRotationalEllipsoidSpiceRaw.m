function [ell_estim,rms]=FindRotationalEllipsoidSpiceRaw(r,fi,lambda,ell)

options=optimset;
options=optimset(options,'UseParallel','always');


[x,y,z]=sph2cart(lambda,fi,r);

rv=[x';y';z'];


[ell_estim,rms]=fminsearch(@(ell) FindParameter2(ell, rv),ell,options);


function SqDist=FindParameter2(ell, r)

[~,H]= cspice_nearpt( r, ell(1), ell(1), ell(2) );

SqDist=sum(sum(H.^2));
ell
