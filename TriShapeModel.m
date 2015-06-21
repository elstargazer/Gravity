
%% Making triangulation of a sphere
[FV,P]=trisphere(6);

x=P(1,:);
y=P(2,:);
z=P(3,:);

[lambda,fi,~]=cart2sph(x,y,z);



%% Load SH shape model
MaxDegreeTopo=360;

load VestaHASTALAVESTAshape_sh720.mat

lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

Degree=lmcosi_shape(:,1);

%% Shape factor
figure; hold on;

xlim([3 max(Degree(:))]);
ylim([0 20]);

gi=ginput(3);

f=interp1(gi(:,1),gi(:,2),Degree,'cubic');

plot(Degree,f,'.');

lmcosi_shape(:,3)=lmcosi_shape(:,3).*f;
lmcosi_shape(:,4)=lmcosi_shape(:,4).*f;

%% Computing radius at verteces
[r,~,~,~]=plm2xyz(lmcosi_shape,fi*180/pi,lambda*180/pi,[],[]);

[xi,yi,zi]=sph2cart(lambda,fi,r');



%% Plotting shape model
figure

tri_fig=trisurf(FV,xi,yi,zi); 
colorbar; axis equal