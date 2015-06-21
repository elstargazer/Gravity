function [fii,lambdai]=MakeGeoGrid(FiStep,LambdaStep)

fi=-90:FiStep:90;
lambda=0:LambdaStep:360;

[fii,lambdai]=meshgrid(fi/180*pi,lambda/180*pi);