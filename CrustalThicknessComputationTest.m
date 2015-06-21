ccc


a1=253.2;
c1=199.3;

a2=280.9;
c2=226.2;


Npoints=10000000;

[fii,lambdai]=GenerateRandomSphericalCoord(Npoints);

x1r=a1.*cos(lambdai).*cos(fii);
y1r=a1.*sin(lambdai).*cos(fii);
z1r=c1.*sin(fii);

x2r=a2.*cos(lambdai).*cos(fii);
y2r=a2.*sin(lambdai).*cos(fii);
z2r=c2.*sin(fii);

r1r=sqrt(x1r.*x1r+y1r.*y1r+z1r.*z1r);
r2r=sqrt(x2r.*x2r+y2r.*y2r+z2r.*z2r);

h1=mean(abs(r2r-r1r));


V1=4/3*pi*a1*a1*c1;
V2=4/3*pi*a2*a2*c2;

R1=(V1*3/4/pi)^(1/3);
R2=(V2*3/4/pi)^(1/3);

4/3*pi*V1^3

h2=abs(R1-R2);


[h1 h2]



(h1-h2)/h2*100