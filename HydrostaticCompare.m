ccc

% Mars
Qv=0.125;
rhoratio=0.486;
eps2=0.00347;

% Neptune
Qv=0.091125;
rhoratio=0.157334;
eps2=0.0254179;

G=6.67e-11;


L=sqrt(2*eps2);

% Neptune
M=102.42e24;
r1=24622000;


r2=r1*(Qv^(1/3));

V=4/3*pi*r1^3;
V2=4/3*pi*r2^2;
V1=V-V2;

rho1=M*rhoratio/(V2+V1*rhoratio);
rho2=M/(V2+V1*rhoratio);

W=sqrt(G*pi*rho1)*L;
T=2*pi/W/3600;

f10=0.1;
f20=0.1;

[fh,fval]=HydrostaticStateExact2l(r1,r2,T,rho1,rho2,f10,f20);

fh'

eh=sqrt(2*fh-fh.*fh);

eh'

