function l=CompensationLevel(R,g,nu,rhomean,rhomantle,rhocrust,E,d,n)


drho=(rhomantle-rhocrust);

D=E.*(d.^3)./(12*(1-nu.^2));

tau=E.*d./((R.^2)*g*drho);

sigma=D./((R.^4).*g*drho);

l=(1+(-3).*(1+2.*n).^(-1).*rhomantle.*rhomean.^(-1)).*((-3).*(1+2.* ...
  n).^(-1).*rhomantle.*rhomean.^(-1)+((-1)+n.*(1+n)+nu).^(-1).*((-1) ...
  +n.*(1+n)+nu+((-4).*n.^2.*(1+n).^2+n.^3.*(1+n).^3).*sigma+((-2)+ ...
  n.*(1+n)).*tau)).^(-1);

