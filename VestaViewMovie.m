close all
ShowVesta

MovieFileName='CeresFlight5.avi';

TimeStart=0;
TimeEnd=3600*5.342*1.4;
TimeStep=10;
Time=TimeStart:TimeStep:TimeEnd;

omega=3.2664e-04;

fi_light=Time.*0;

lambda_light=Time*0;

[x_light,y_light,z_light]=sph2cart(lambda_light,fi_light,1);

% 
% Rorbit=460000;

mu=17.28e9;

% n=sqrt(mu/(Rorbit^3));
% 
% x_sc=Rorbit*cos(n*Time)+rand(1)/100;
% z_sc=Rorbit*sin(n*Time)+rand(1)/100;
% y_sc=0.*x+rand(1)/100;

a=1000000;
e=0.75;
i=85/180*pi;
w=270/180*pi;
W=225/180*pi;
M0=180;
t0=0;
t=Time;

[x_sc,y_sc,z_sc,vx_sc,vy_sc,vz_sc]=el_koor(a,e,i,w,W,M0,t0,t,mu);

i=1;

campos([x_sc(i),y_sc(i),z_sc(i)]);

set(gca,'CameraTarget',[0 0 0],'CameraViewAngle',65)
set(gca,'CameraUpVectorMode','manual','CameraPositionMode','manual',...
    'CameraTargetMode','manual','CameraViewAngleMode','manual');

NumberOfFrames=numel(Time);

vidObj = VideoWriter(MovieFileName);
open(vidObj);
axis vis3d

% ArrowLength=500000;
% line([0 ArrowLength],[0 0],[0 0],'color','r');
% line([0 0],[0 ArrowLength],[0 0],'color','g');
% line([0 0],[0 0],[0 ArrowLength],'color','b');

for i=1:NumberOfFrames
   
set(light_handle,'Position',rot(omega*Time(i),3)*[x_light(i) y_light(i) z_light(i)]');


r_star_new=rot(omega*Time(i),3)*r_star';
set(plot_star,'XData',r_star_new(1,:));
set(plot_star,'YData',r_star_new(2,:));
set(plot_star,'ZData',r_star_new(3,:));


r_sc=[x_sc(i);y_sc(i);z_sc(i)];
r_cam=rot(omega*Time(i),3)*r_sc;
campos(r_cam');


rv_sc_veldir=[vx_sc(i);vy_sc(i);vz_sc(i)];
rv_sc_veldir=rv_sc_veldir/norm(rv_sc_veldir);
rv_cam_up=rot(omega*Time(i),3)*rv_sc_veldir;
% set(gca,'CameraUpVector',rv_cam_up')

camup(rv_cam_up);

%  drawnow

currFrame = getframe(gcf);

writeVideo(vidObj,currFrame);

i/NumberOfFrames
end

close(vidObj);

% rv_sc=[x_sc(i);y_sc(i);z_sc(i)];
% 
% rv_sc_veldir=[y_sc(i);-x_sc(i);0];
% rv_sc_veldir=rv_sc_veldir/norm(rv_sc_veldir);
% 
% rv_cam=rot(omega*Time(i),3)*rv_sc;
% rv_cam_up=rot(omega*Time(i),3)*rv_sc_veldir;
% 
% campos(r_cam');

% set(gca,'CameraUpVector',rv_cam_up');