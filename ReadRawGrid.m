function [ri,fii,lambdai]=ReadRawGrid(FileName,Var)

finfo = ncinfo(FileName);

lambda = ncread(FileName,Var{1});
fi = ncread(FileName,Var{2});
ri = ncread(FileName,Var{3});

[fii,lambdai]=meshgrid(fi/180*pi,lambda/180*pi);