function mu=muEllipsoid(a,b,c,rho)

V=4/3*pi*a*b*c;
m=rho*V;
G=6.67384e-11;
mu=G*m;