function Z=AiryAdmittanceFun(x,n)

g=1.6;
nu=0.25;

% rhocrust=2600;
rhomantle=3220;
rhomean=3344;
Rref=1738000;

% h_comp=2*x(2);
% (R,g,nu,rhomean,rhomantle,rhocrust,E,d,n,tl)
Z=AiryAdmittance(Rref,g,nu,rhomean,rhomantle,x(3),x(1),x(2),n,2*x(2));

