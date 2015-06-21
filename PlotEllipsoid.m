function PlotEllipsoid(ag,bg,cg)

FiStep=10;
LambdaStep=10;
[fii,lambdai]=meshgrid( -pi/2:FiStep/180*pi:pi/2 , 0:LambdaStep/180*pi:2*pi);
x=ag*cos(lambdai).*cos(fii);
y=bg*sin(lambdai).*cos(fii);
z=cg*sin(fii);

surf(x,y,z);
