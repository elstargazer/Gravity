ccc

% VestaMatFileName='VestaHASTALAVESTAshape_sh1000.mat';
VestaMatFileName='~/Dawn/Balmino2/VestaTest/SH_VestaHASTALAVESTAshape_6min';

MarsFileName='/Users/antonermakov/Mars/MarsTopo719.shape';


%% Vesta topo sh

% load(VestaMatFileName);lmcosi=lmcosi_shape;
load VestaHASTALAVESTAshape_sh1500.mat

lmcosi_Vesta=lmcosi_shape;
clear lmcosi_shape
%% Planets topo sh

lmcosi_Mars=load(MarsFileName);


MaxDegreeTopo=400;
lmcosi_Mars=TruncateGravityModel(lmcosi_Mars,MaxDegreeTopo,1);
lmcosi_Vesta=TruncateGravityModel(lmcosi_Vesta,MaxDegreeTopo,1);

lmcosi=lmcosi_Vesta;

%% Create a mesh
Step=.1;

[lambdai,fii]=meshgrid(0:Step:360,90:-Step:-90);

lambdai=lambdai/180*pi;
fii=fii/180*pi;



ri1=plm2xyz(lmcosi_Vesta,Step);
[xi1,yi1,zi1]=sph2cart(lambdai,fii,ri1);

ri2=plm2xyz(lmcosi_Mars,Step);
[xi2,yi2,zi2]=sph2cart(lambdai,fii,ri2);


figure('Position',[1 1 1500 1500]);


xi=xi1;
yi=yi1;
zi=zi1;
ri=ri1;

plot1=surf(xi,yi,zi,ri);

lighting phong
shading interp
axis tight equal off
camlight('headlight');

view(20,-20);

% Prepare the new file.
vidObj = VideoWriter('VestaToMars.avi');
open(vidObj);


tic



for k=0:(1/500):1
    
%     lmcosi(:,3)=k*lmcosi_Mars(:,3)+(1-k)*lmcosi_Vesta(:,3);
%     lmcosi(:,4)=k*lmcosi_Mars(:,4)+(1-k)*lmcosi_Vesta(:,4);
%     
%     
%     
%     ri=plm2xyz(lmcosi,Step);
%     
%     
%     [xi,yi,zi]=sph2cart(lambdai,fii,ri);

    xi=(1-k)*xi1+k*xi2;
    yi=(1-k)*yi1+k*yi2;
    zi=(1-k)*zi1+k*zi2;
    ri=(1-k)*ri1+k*ri2;
    
    
    
    set(plot1,'XData',xi);
    set(plot1,'YData',yi);
    set(plot1,'ZData',zi);
    set(plot1,'CData',ri);
    
    l=2900000;
    axis([min(xi(:)) max(xi(:)) min(yi(:)) max(yi(:)) min(zi(:)) max(zi(:))  ]);
    
   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
    
   k
    
end

% Close the file.
close(vidObj);

toc


ccc

SlopeMapTest