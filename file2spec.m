function [sdl,l] = file2spec(filename,L)

data = load(filename);

x = data(:,1);
y = data(:,2);
z = data(:,3);

[lon,lat,r] = cart2sph(x,y,z);

lmcosi = xyz2plm(r,L,'irr',lat*180/pi,lon*180/pi);

[sdl,l] = plm2spec(lmcosi);






