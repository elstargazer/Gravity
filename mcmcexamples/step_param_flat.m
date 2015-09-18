function param_new=step_param_flat(param,s)

r1=param(1);
rho1=param(2);

% rho2=-(3*M-4*pi*(r1.^3).*rho1)./(4*pi*(r1.^3)-4*pi*(r2^3));

r1=r1+randn*s(1);
rho1=rho1+randn*s(2);

param_new=[r1 rho1];
