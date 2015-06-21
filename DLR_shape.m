ccc

filename='/Users/antonermakov/Dawn/Vesta_Shape_Models/DLR/2014/global_64.cub';

[label, core, isis_version] = read_isis_cub(filename);

core(1,:)=[];
core(end,:)=[];
core(:,1)=[];
core(:,end)=[];


step=1/64;

s=size(core);

% lambda=(0:step:360)/180*pi;
% fi=(90:-step:-90)/180*pi;

lambda=linspace(0,360,s(2));
fi=linspace(-90,90,s(1));

[lambdai,fii]=meshgrid(lambda,fi);

[x,y,z]=sph2cart(lambdai/180*pi,fii/180*pi,core);


% surf(x(1:10:end,1:10:end),y(1:10:end,1:10:end),z(1:10:end,1:10:end),core(1:10:end,1:10:end)); shading interp;
% StandardLight

core(core<0)=NaN;

% clear x y z lambda fi

% WriteXYZ(lambdai,fii,core,'VestaShapeDLR64.xyz');

inr=1:300;

[xs,ys]=sph2stereo(lambdai(inr,:)/180*pi,fii(inr,:)/180*pi,0,-pi/2,265000);


figure; hold on;
surf(xs,ys,core(inr,:)); shading interp;
colorbar


% plot(xs,ys,'.')


figure; hold on;
surf(x(inr,1:1:end),y(inr,1:1:end),z(inr,1:1:end)); shading interp;
colorbar


% 
% 
% 
% figure('Position',[1 1 1000 500]); hold on;
% set(gca,'FontSize',20);
% 
% % Prepare the new file.
% vidObj = VideoWriter('VestaDLRSection.avi');
% open(vidObj);
% 
% ps=plot(circshift(lambda,[0 fix(s(2)/2)]),circshift(core(1,:),[0 fix(s(2)/2)]),'.');
% 
% xlim([0 360]);
% set(gca,'XTick',0:30:360);
% 
% currFrame = getframe(gcf);
% writeVideo(vidObj,currFrame);
% 
% xlabel('Longitude [deg]','FontSize',20);
% ylabel('Elevation [m]','FontSize',20);
% title(['Latitude = ' num2str(fi(1))  ' deg'],'FontSize',20);
% 
% for i=2:s(1)
% 
% set(ps,'YData',circshift(core(i,:),[0 fix(s(2)/2)]));
% title(['Latitude = ' num2str(fi(i))  ' deg'],'FontSize',20);
% 
% currFrame = getframe(gcf);
% writeVideo(vidObj,currFrame);
% 
% 
% end
%  
%  
% close(vidObj);