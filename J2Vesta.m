function J2=J2Vesta(r1,r2,f2,rho1,rho2,J2outer,Rref,M)

a2=r2/((1-f2).^(1/3));
c2=r2/((1-f2).^(1/3))-f2*r2/((1-f2).^(1/3));

J2=((c2*c2-a2*a2)*(rho2-rho1)*(4/3*pi*r2^3)/(5*Rref*Rref)-J2outer*sqrt(5)*rho1*(4/3*pi*r1^3))/(-sqrt(5)*M);