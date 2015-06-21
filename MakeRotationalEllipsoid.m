function [x,y,z]=MakeRotationalEllipsoid(ag,bg,FiStep,LambdaStep)

[lambdai,fii]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2);
x=ag*cos(lambdai).*cos(fii);
y=ag*sin(lambdai).*cos(fii);
z=bg*sin(fii);