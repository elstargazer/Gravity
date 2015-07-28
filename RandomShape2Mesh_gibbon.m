ccc

L    = 100;
beta = -3;

[lmcosi,bta,bto,sdl,el]=plm2rnd(L,beta,1,3);

[sdl,l] = plm2spec(lmcosi);

figure; hold on;
set(gca,'XScale','log');
set(gca,'YScale','log');
plot(l,sdl,'-r');

[r] = plm2xyz(lmcosi,lat,lon);





