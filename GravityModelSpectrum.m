function GravityModelSpectrum(lmcosi,color)

[sdl,l,~,~,~,~]=plm2spec(lmcosi);
[sdl_err,l_err,~,~,~,~]=plm2spec([lmcosi(:,1:2) lmcosi(:,5:6)]);


set(gca,'FontSize',20,'YScale','log');


plot(l,sdl,'o-','LineWidth',2,'Color',color);
plot(l_err,sdl_err,'o--','LineWidth',2,'Color',color);

xlabel('Degree','FontSize',20);
ylabel('Coef. power','FontSize',20);