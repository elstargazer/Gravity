function Z=AiryAdmittance(R,g,nu,rhomean,rhomantle,rhocrust,E,d,n,tl)

l=CompensationLevel(R,g,nu,rhomean,rhomantle,rhocrust,E,d,n);

% l=1;

Coef=1.6*1e5/R*1000;

Z=Coef*3./(2.*n+1).*(rhocrust./rhomean).*(1-((1-tl./R).^n).*l).*(n+1);