ccc

minCoef=0.992;
maxCoef=1.008;
N=10;
Coef=linspace(minCoef,maxCoef,N);

progressbar(0)



for i=1:numel(Coef)

[~,~,~,rms(i)]=GravityPotentialTaylorTestOptim(Coef(i));

progressbar(i/numel(Coef))


end

progressbar(1);

plot(Coef,rms);

P=polyfit(Coef,rms,2);

DP = polyder(P);


CoefOpt=roots(DP);

[fii,lambdai,ri_g,rms]=GravityPotentialTaylorTestOptim(0.997872488650079);

