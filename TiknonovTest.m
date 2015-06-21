close all



vidObj = VideoWriter('Lambda.avi');
open(vidObj);


l=10.^(2:-0.25:-8);

MaxDegree=15;
Resolution=1;



AGUaxes
cbar=colorbar('FontSize',20);
ylabel(cbar,'Density [kg/m^{3}] ','FontSize',20);

for i=1:numel(l)

lambda=l(i);


Astar=(A+lambda*(G'*G));


drho1=Astar\b;


for n=1:Nelems/kp
rho1_s(n:Nelems/kp:Nelems)=rho1(n:Nelems/kp:Nelems)+drho1(n);
end


lmcosi_den=xyz2plm(rho1_s,MaxDegree,'irr',fi1m'*180/pi,lambda1m'*180/pi);
[rho1_s_sh,lon_sh,lat_sh]=plm2xyz(lmcosi_den,Resolution);

[rho1_s_shm,lon_shm,lat_shm]=plm2xyz(lmcosi_den,fi1m*180/pi,lambda1m*180/pi);

[lon_sh,lat_sh]=meshgrid(lon_sh/180*pi,lat_sh/180*pi);


pcolorm(lat_sh,lon_sh,rho1_s_sh); shading interp;
% contourm(lat_sh,lon_sh,rho1_s_sh,'LevelStep',100,'ShowText','on','Color','k');


title(['\lambda = ' num2str(lambda)],'FontSize',20);
caxis([2300 3600]);


currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);


unplot

i
lambda


end

close(vidObj);