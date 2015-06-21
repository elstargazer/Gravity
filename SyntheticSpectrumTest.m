ccc

%% Creating artificial field
MaxDegree=40;
BTA=-1;

[lmcosi,bta,bto,sdl,el]=plm2rnd(MaxDegree,BTA,2);


S_global=CrossPower(lmcosi,lmcosi);


r=plm2xyz(lmcosi,1);

% MapRadialGrid(r);

%% Localization
L=5;


circle_rad=50;
lon_center=100;
lat_center=45;


[G,V,N,J]=glmalphapto(circle_rad,L,lon_center,90-lat_center);

s=size(G);


for i=1:s(2)       
    lmcosi_window{i}=glm2lmcosi(G,i);
    [rl(:,:,i),lor,lar,Plm]=plm2xyz(lmcosi_window{i},1);
    [lmcosi_w{i},~]=xyz2plm(rl(:,:,i).*r*3,MaxDegree,'im');
    S_local(:,i)=CrossPower(lmcosi_w{i},lmcosi_w{i});
end

%% plot eigenvalue

figure
plot(V);


Tapers=find(V>0.9);
%% plot power spectra

figure; hold on;

semilogy(S_global,'-b','LineWidth',3);


S_local_mean=mean(S_local(:,Tapers),2);


S_local_std=std(S_local(:,Tapers),0,2);

semilogy(S_local(:,Tapers),'--r','LineWidth',1);
semilogy(S_local_mean,'-r','LineWidth',3);
semilogy(S_local_mean+S_local_std,'--r','LineWidth',2);
semilogy(S_local_mean-S_local_std,'--r','LineWidth',2);


