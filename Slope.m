function [slope,slope_az]=Slope(xi,yi,zi)

% nxi=zeros(size(xi));
% nyi=zeros(size(yi));
% nzi=zeros(size(zi));

%for i=1:size(nxi,2)
%    [nxi(:,i),nyi(:,i),nzi(:,i)] = surfnorm(xi,yi,zi);
%end

[nxi,nyi,nzi] = surfnorm(xi,yi,zi);

% 
% figure; hold on;
% surf(xi,yi,zi);
% quiver3(xi(1:16:end,1:16:end),yi(1:16:end,1:16:end),zi(1:16:end,1:16:end),...
%     nxi(1:16:end,1:16:end),nyi(1:16:end,1:16:end),nzi(1:16:end,1:16:end),10);
% 
% axis equal;

nxi=-nxi;
nyi=-nyi;
nzi=-nzi;

mean_nxi=(nxi(:,1)+nxi(:,end))/2;
mean_nyi=(nyi(:,1)+nyi(:,end))/2;
mean_nzi=(nzi(:,1)+nzi(:,end))/2;

mean_norm=sqrt(mean_nxi.^2+mean_nyi.^2+mean_nzi.^2);

mean_nxi=mean_nxi./mean_norm;
mean_nyi=mean_nyi./mean_norm;
mean_nzi=mean_nzi./mean_norm;

nxi(:,1)=mean_nxi;
nyi(:,1)=mean_nyi;
nzi(:,1)=mean_nzi;

nxi(:,end)=mean_nxi;
nyi(:,end)=mean_nyi;
nzi(:,end)=mean_nzi;

mean_nxi=mean(nxi(1,:));
mean_nyi=mean(nyi(1,:));
mean_nzi=mean(nzi(1,:));

nxi(1,:)=nxi(1,:)*0+mean_nxi;
nyi(1,:)=nyi(1,:)*0+mean_nyi;
nzi(1,:)=nzi(1,:)*0+mean_nzi;

mean_nxi=mean(nxi(end,:));
mean_nyi=mean(nyi(end,:));
mean_nzi=mean(nzi(end,:));

mean_norm=sqrt(mean_nxi.^2+mean_nyi.^2+mean_nzi.^2);

mean_nxi=mean_nxi/mean_norm;
mean_nyi=mean_nyi/mean_norm;
mean_nzi=mean_nzi/mean_norm;

nxi(end,:)=nxi(end,:)*0+mean_nxi;
nyi(end,:)=nyi(end,:)*0+mean_nyi;
nzi(end,:)=nzi(end,:)*0+mean_nzi;

ag=1;
bg=1;

[s_up,s_east,s_north]=GravityComponents(nxi,nyi,nzi,xi,yi,zi,ag,bg);

slope_az=atan2(s_north,s_east);
slope=acos(s_up)*180/pi;

slope=real(slope);






