function lmcosi_grad=ComputeGradientCoefficients(lmcosi,DerOrder,r)

lmcosi_grad=lmcosi;

Coef=RadialGradientFactor(lmcosi(:,1),DerOrder,r);

lmcosi_grad(:,3)=lmcosi(:,3).*Coef;
lmcosi_grad(:,4)=lmcosi(:,4).*Coef;