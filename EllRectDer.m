function [d11,d12,d13,d21,d22,d23,d31,d32,d33]=EllRectDer(lambda_ell,x,y,z,a,b,c)

x=x(:);
y=y(:);
z=z(:);

h2 = (a^2-b^2);
k2 = (a^2-c^2);

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

lambda=lambda_ell(:,1);
mu=lambda_ell(:,2);
nu=lambda_ell(:,3);

lambda2=lambda.*lambda;
mu2=mu.*mu;
nu2=nu.*nu;

d11=x.*(lambda2-k2).*(lambda2-h2)./(lambda.*(lambda2-mu2).*(lambda2-nu2));
d12=x.*(mu2-k2).*(mu2-h2)./(mu.*(mu2-lambda2).*(mu2-nu2));
d21=y.*lambda.*(lambda2-k2)./((lambda2-mu2).*(lambda2-nu2));
d22=y.*mu.*(mu2-k2)./((mu2-lambda2).*(mu2-nu2));
d23=y.*nu.*(nu2-k2)./((nu2-lambda2).*(nu2-mu2));
d31=z.*lambda.*(lambda2-h2)./((lambda2-mu2).*(lambda2-nu2));
d32=z.*mu.*(mu2-h2)./((mu2-lambda2).*(mu2-nu2));
d33=z.*nu.*(nu2-h2)./((nu2-lambda2).*(nu2-mu2));

Cond=(nu==0);
d13(Cond)=h.*k./(lambda(Cond).*mu(Cond));
d13(~Cond)=x(~Cond).*(nu2(~Cond)-k2).*(nu2(~Cond)-h2)./...
    (nu(~Cond).*(nu2(~Cond)-lambda2(~Cond)).*(nu2(~Cond)-mu2(~Cond)));

d13=d13';