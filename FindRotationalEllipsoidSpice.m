 function [ell_estim,rms]=FindRotationalEllipsoidSpice(ri,fii,lambdai,ell)

 Npoints=100000;

x_rand=rand(1,Npoints)-0.5;
y_rand=rand(1,Npoints)-0.5;
z_rand=rand(1,Npoints)-0.5;

r_rand=sqrt(x_rand.^2+y_rand.^2+z_rand.^2);

x_rand=x_rand(r_rand<0.5);
y_rand=y_rand(r_rand<0.5);
z_rand=z_rand(r_rand<0.5);


[lambda_rand,fi_rand,~]=cart2sph(x_rand,y_rand,z_rand);
lambda_rand=lambda_rand+2*pi.*(lambda_rand<0);

r_rand=griddata(lambdai,fii,ri,lambda_rand,fi_rand,'nearest'); 

options=optimset;
options=optimset(options,'UseParallel','always');

[x,y,z]=sph2cart(lambda_rand,fi_rand,r_rand);

x=x(:);
y=y(:);
z=z(:);

r=[x,y,z]';

[ell_estim,rms]=fminsearch(@(ell) FindParameter2(ell,r),ell,options);


function SqDist=FindParameter2(ell, r)

[~,H]= cspice_nearpt( r, ell(1), ell(1), ell(2) );

SqDist=sum(sum(H.^2));
ell
