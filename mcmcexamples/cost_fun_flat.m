function cost=cost_fun_flat(param,data,T,r2,M)

r1=param(1);
rho1=param(2);

rho2=-(3*M-4*pi*(r1.^3).*rho1)./(4*pi*(r1.^3)-4*pi*(r2^3));

% f=HydrostaticStateExact2l(r2,r1,T,rho2,rho1,0.1,0.1);
f=HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,0.1,0.1,2);

f_obs=data(1);
J2_obs=data(2);

J2=RadFlat2J2(r2,r1,f(2),f(1),rho2,rho1,476000);

cost=((f(2) - f_obs).^2)./(0.004^2) + ...
     ((J2  - J2_obs).^2)./(25e-5^2);



