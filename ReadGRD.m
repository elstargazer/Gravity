function [lambdai,fii,ri]=ReadGRD(filename)

finfo = ncinfo(filename);

lambdai=ncread(filename,finfo.Variables(1).Name);
fii=ncread(filename,finfo.Variables(2).Name);
ri=ncread(filename,finfo.Variables(3).Name)';

% lambdai=ncread(filename,'x');
% fii=ncread(filename,'y');
% ri=ncread(filename,'z')';

[lambdai,fii]=meshgrid(lambdai,fii);