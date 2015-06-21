% ro=3.458942584995802e+03;
% G=6.67384e-11;
% mu=17.2883080293e9;
% Rref=265000;

ro=3346.4;
G=6.67384e-11;
mu=4902.7779e9;
Rref=1737000;

vidObj = VideoWriter('TopoPowerCorrelation-Moon.avi');
open(vidObj);


MoonFileName='/Users/antonermakov/GRAIL/Topography/SH/LRO_LTM05_2050_SHA.TAB.txt';

lmcosi_full=load(MoonFileName);

% plotplm(lmcosi,[],[],2,1,[],[],[]);


for MaxTopoPower=1:20

MaxDegreeTopo=300;
MaxDegreeGrav=300;
Resolution=1;

lmcosi=TruncateGravityModel(lmcosi_full,MaxDegreeTopo,1);


% lmcosi(4,3)=0;

[ri_shape,~]=plm2xyz(lmcosi,Resolution);

lmcosi_gt=TopoSH2GravitySH(ri_shape,mu,ro,Rref,MaxDegreeTopo,MaxDegreeGrav,MaxTopoPower);

lmcosi_gt(1,:)=[];
lmcosi(1,:)=[];


SphericalHarmonicCorrelation(lmcosi,lmcosi_gt,'m','-');


   axis([2 MaxDegreeTopo 0 1]);
   
   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
   unplot
   
end   
   
close(vidObj);