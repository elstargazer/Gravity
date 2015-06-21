ccc
filename='/Users/ermakov/Dawn/ShapeModel/Gaskell/20140513/MAPINFO.TXT';


in=fopen(filename,'r');
data=textscan(in,'%s %f %f %f %f %f %f %f %f %f%*[^\n]');

lat=data{3};
lon=data{4};
rad=data{5};

std=data{9};

AGUaxes;

scatterm(lat(1:10:end)/180*pi,lon(1:10:end)/180*pi,50,std(1:10:end),'filled');

colorbar

caxis([0 50]);