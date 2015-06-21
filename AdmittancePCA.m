%% PCA of Admittance
% [coeff,score,latent,tsquared,explained,mu] = pca(Z_t);

 %Z_t=density_t;

 sz=size(Z_t);

[w,pc,ev,tsquared] = princomp(Z_t);

mean_Z_t=mean(Z_t);
mean_Z_t=repmat(mean_Z_t,[sz(1) 1]);

Ncomp=3;

RGBcomp=[1 2 3];
MaxDegreeExp=24;
% Npc=5;
for Npc=1:Ncomp
% Y=w(:,Npc);
Y=w(:,Npc)';

lmcosi_Y=xyz2plm(Y,MaxDegreeExp,'irr',fii,lambdai);
[Y_sh(:,:,Npc),lonsh,latsh]=plm2xyz(lmcosi_Y,Resolution);
[lonsh,latsh]=meshgrid(lonsh,latsh);

AGUaxes
hp=pcolorm(latsh/180*pi,lonsh/180*pi,Y_sh(:,:,Npc)); shading interp
cb1=colorbar('FontSize',20);
% ylabel(cb1,'Effective density [kg/m ^3] ','FontSize',20);
 %PlotVestaFeatures
title(['Principal component # ' num2str(Npc)],'FontSize',20);

end

s_sh=size(Y_sh);
RGB=zeros(s_sh(1),s_sh(2),3);
RGB(:,:,1)=Y_sh(:,:,RGBcomp(1));
RGB(:,:,2)=Y_sh(:,:,RGBcomp(2));
RGB(:,:,3)=Y_sh(:,:,RGBcomp(3));

Rrange= max(max(RGB(:,:,1)))-min(min(RGB(:,:,1)));
Grange= max(max(RGB(:,:,2)))-min(min(RGB(:,:,2)));
Brange= max(max(RGB(:,:,3)))-min(min(RGB(:,:,3)));

RGB(:,:,1)=(RGB(:,:,1)-min(min(RGB(:,:,1))))/Rrange;
RGB(:,:,2)=(RGB(:,:,2)-min(min(RGB(:,:,2))))/Grange;
RGB(:,:,3)=(RGB(:,:,3)-min(min(RGB(:,:,3))))/Brange;


 AGUaxes
% axesm('mollweid','frame','on','FontSize',24);
pl_pcr=surfm(latsh/180*pi,lonsh/180*pi,ri); shading interp
% colormap jet

% plot_pc=surflsrm(latsh/180*pi,lonsh/180*pi,ri);
%  
set(pl_pcr,'CData',RGB);
%PlotVestaFeatures

ShowVestaLow
lighting phong

light_handle1=light('Style','infinite');
light_handle2=light('Style','infinite');


set(h,'CData',RGB);
lambda_light=0;
fi_light=90;
[x_light, y_light, z_light]=sph2cart(lambda_light/180*pi,fi_light/180*pi,1);
set(light_handle1,'Position',[x_light, y_light, z_light]);

lambda_light=0;
fi_light=-90;
[x_light, y_light, z_light]=sph2cart(lambda_light/180*pi,fi_light/180*pi,1);
set(light_handle2,'Position',[x_light, y_light, z_light]);
%% plot principal componet
take_comp=[1 2 3 4 5];

figure; hold on;
set(gca,'FontSize',20);
plot(GoodDegrees,pc(:,take_comp));
% legend(take_comp)

%% Reconstruct data
take_comp=[1 2 3 4 5 6];

Z_tr=pc(:,take_comp)*w(:,take_comp)';
Z_tr=Z_tr+mean_Z_t;

for i=1:sz(2)
%     density_mean_tr(i)=mean((Z_t(:,i))./Z_t_shape(:,i)*rho);
%     density_mean_tr(i)=mean(Z_tr(:,i)./Z_t_shape(:,i)*rho);
%     p_tr(i,:)=polyfit(GoodDegrees,(Z_tr(:,i)./Z_t_shape(:,i)*rho)',1);
     density_mean_tr(i)=mean(Z_tr(:,i));
     p_tr(i,:)=polyfit(GoodDegrees,(Z_tr(:,i))',1);
%      density_std(i)=std(Z_tr(:,i));

end

%
lmcosi_den_tr=xyz2plm(density_mean_tr,MaxDegreeExp,'irr',fii,lambdai);
[density_mean_sh_tr,lonsh,latsh]=plm2xyz(lmcosi_den_tr,Resolution);
[lonsh,latsh]=meshgrid(lonsh,latsh);

AGUaxes
pcolorm(latsh/180*pi,lonsh/180*pi,density_mean_sh_tr); shading interp
cb1=colorbar('FontSize',20);
ylabel(cb1,'Effective density [kg/m ^3] ','FontSize',20);
PlotVestaFeatures
title('Effective density ','FontSize',20);
%
lmcosi_p1_tr=xyz2plm(p_tr(:,1),MaxDegreeExp,'irr',fii,lambdai);
[p1_sh_tr,lonsh,latsh]=plm2xyz(lmcosi_p1_tr,Resolution);
[lonsh,latsh]=meshgrid(lonsh,latsh);

AGUaxes
pcolorm(latsh/180*pi,lonsh/180*pi,p1_sh_tr); shading interp
cb1=colorbar('FontSize',20);
ylabel(cb1,'Effective density slope [kg/m ^3/degree] ','FontSize',20);
PlotVestaFeatures
title('Effective density slope ','FontSize',20);

%
% lmcosi_std_tr=xyz2plm(density_std,MaxDegreeExp,'irr',fii,lambdai);
% [density_std_sh_tr,lonsh,latsh]=plm2xyz(lmcosi_std_tr,Resolution);
% [lonsh,latsh]=meshgrid(lonsh,latsh);
% 
% AGUaxes
% pcolorm(latsh/180*pi,lonsh/180*pi,density_std_sh_tr); shading interp
% cb1=colorbar('FontSize',20);
% ylabel(cb1,'Effective std [kg/m ^3] ','FontSize',20);
% PlotVestaFeatures
% title('Effective std ','FontSize',20);

100*ev(take_comp)/sum(ev)

figure

bar(100*ev(take_comp)/sum(ev))


