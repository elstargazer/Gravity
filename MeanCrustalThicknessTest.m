
ccc

tic

ro_core=7400;

d_ro=50;

ro_crust=2700:d_ro:3100;
ro_mantle=3000:d_ro:3500;

[ro_crusti,ro_mantlei]=meshgrid(ro_crust,ro_mantle);


Condition=ro_mantlei<ro_crusti;
ro_crusti(Condition)=NaN;
ro_mantlei(Condition)=NaN;

s=size(ro_crusti);

progressbar(0)




[C_gt,S_gt]=ReadBalminoSH('file_harmo_pot');
lmcosi_gt=CS2lmcosi(C_gt,S_gt);

Rref=265000;
% C_gt(1,2)=offset(1)/Rref;
% S_gt(1,2)=offset(2)/Rref;
% C_gt(1,1)=offset(3)/Rref;

GravityFileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.SHA';
[lmcosi_obs_full,Rref,mu,mu_std]=ReadGRAILGravityModel(GravityFileName);

[C_obs,S_obs,C_obs_std,S_obs_std,mu_obs,Rref]=LoadGravityModel(GravityFileName);



for i=1:s(1)
    for j=1:s(2)
        if (ro_crusti(i,j)~=NaN)
            
           [mean_crustal_thickness(i,j),max_crust(i,j),min_crust(i,j),max_bouguer(i,j),min_bouguer(i,j),f_mantle_new(i,j)]...
               =MeanCrustalThickness(ro_crusti(i,j),ro_mantlei(i,j),ro_core,C_obs,S_obs,mu_obs,Rref,lmcosi_obs_full,C_gt,S_gt,lmcosi_gt);
           
           clc
           [ro_crusti(i,j) ro_mantlei(i,j)]
        end
        
        progressbar(i/s(1));
        
        fclose all;
    end
end

bouguer_range=max_bouguer-min_bouguer;

%% Bouguer range

% figure
% h0=pcolor(ro_crusti,ro_mantlei,bouguer_range);
% colorbar
% shading interp
% hold on;
% alpha(0.5)
% axis([2700 3050 3150 3500]);

%%


%% interpolation

ro_crusti=ro_crusti(:);
ro_mantlei=ro_mantlei(:);
mean_crustal_thickness=mean_crustal_thickness(:);

min_crust=min_crust(:);
max_crust=max_crust(:);
bouguer_range=bouguer_range(:);


Condition2=isnan(ro_crusti);

ro_crusti(Condition2)=[];
ro_mantlei(Condition2)=[];
mean_crustal_thickness(Condition2)=[];

min_crust(Condition2)=[];
max_crust(Condition2)=[];

bouguer_range(Condition2)=[];


d_ro2=2;

[ro_crustii,ro_mantleii]=meshgrid(min(ro_crust(:)):d_ro2:max(ro_crust(:)),min(ro_mantle(:)):d_ro2:max(ro_mantle(:)));


% min_crust(min_crust<-200)=NaN;

mean_crustal_thickness(mean_crustal_thickness>150)=NaN;



mean_crustal_thicknessi=griddata(ro_crusti,ro_mantlei,mean_crustal_thickness,ro_crustii,ro_mantleii,'cubic');

min_crusti=griddata(ro_crusti,ro_mantlei,min_crust,ro_crustii,ro_mantleii,'cubic');
max_crusti=griddata(ro_crusti,ro_mantlei,max_crust,ro_crustii,ro_mantleii,'cubic');


MaxBouguerLimit=3000;

bouguer_range(bouguer_range>MaxBouguerLimit)=NaN;


bouguer_rangei=griddata(ro_crusti,ro_mantlei,bouguer_range,ro_crustii,ro_mantleii,'cubic');


mean_crustal_thickness(mean_crustal_thickness<0)=NaN;

%% Making a figure

figure


v=0:10:100;

[C,h1]=contour(ro_crustii,ro_mantleii,mean_crustal_thicknessi,v,'LineWidth',2,'Color','k');

clabel(C,h1,v);

xlabel('\rho_{crust}','FontSize',12);
ylabel('\rho_{mantle}','FontSize',12);

% set(gca, 'Position',[0 0 1 1])
set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

hold on;

% print(gcf, '-dpsc', 'ContoursCrustalThickness_ro_core7800.eps');

toc

%% Max thickness
% figure
% v=-100:10:200;
% 
% [C,h2]=contour(ro_crusti,ro_mantlei,max_crust,v,'LineWidth',2,'Color','k');
% 
% clabel(C,h2,v);
% 
% 
% xlabel('\rho_{crust}','FontSize',12);
% ylabel('\rho_{mantle}','FontSize',12);


%% Min thickness
% figure
% v=-100:10:100;
% 
% [C,h3]=contour(ro_crusti,ro_mantlei,min_crust,v,'LineWidth',2,'Color','k');
% 
% clabel(C,h3,v);
% 
% 
% xlabel('\rho_{crust}','FontSize',12);
% ylabel('\rho_{mantle}','FontSize',12);

%%
contour(ro_crustii,ro_mantleii,min_crusti,0*[1 1],'LineWidth',2,'Color','r','LineStyle','--');


caxis([0 max(min_crust)])



% pcolor(ro_crustii,ro_mantleii,min_crusti);shading interp

contour(ro_crustii,ro_mantleii,max_crusti,99.9:100.1,'LineWidth',2,'Color','b');

v4=100:10:800;

% [C4,h4]=contour(ro_crustii,ro_mantleii,bouguer_rangei,v4,'LineWidth',2,'Color','g');

% clabel(C4,h4,v4);


axis([2700 3000 3150 3500]);


% save('CrustalThickessDataSet_ro_core7100.mat');

PrintWhite('CrustalThickessDataSet_ro_core7400_final.eps')



