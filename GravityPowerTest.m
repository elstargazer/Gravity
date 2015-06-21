%% Load gravity models
filename1='../GravityModel/1200/JGGRAIL_1200C12A_SHA.TAB';
filename2='../GravityModel/1200/JGGRAIL_1200C12B_SHA.TAB';


[lmcosi_grav_a,Rref,mu,mu_std]=ReadGRAILGravityModel(filename1);
[lmcosi_grav_b,Rref,mu,mu_std]=ReadGRAILGravityModel(filename2);

lmcosi_grav_diff=lmcosi_grav_a;
lmcosi_grav_diff(:,3:4)=lmcosi_grav_a(:,3:4)-lmcosi_grav_b(:,3:4);


%% Compute power spectra
[l,sdl_a]=plm2spec(lmcosi_grav_a);
[l,sdl_b]=plm2spec(lmcosi_grav_b);
[l,sdl_diff]=plm2spec(lmcosi_grav_diff);

figure; hold on;

plot(l,sdl_a,'r-');
plot(l,sdl_b,'b-');
plot(l,sdl_diff,'b-');

xlabel('Degree []','FontSize',12);
ylabel('Spectral power []','FontSize',12);
set(gca,'YScale','log');