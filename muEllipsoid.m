function mu=muEllipsoid(a,b,c,ro)

V=4/3*pi*a*b*c;
m=ro*V;
G=6.67384e-11;
mu=G*m;