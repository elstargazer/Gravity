function [x,y,z]=ReadSPCGrid(filename,step)

Npts=10000;
[x_rand,y_rand,z_rand]=LoadOpNavShape(filename,Npts,'rand');
[x_grid,y_grid,z_grid]=LoadOpNavShape(filename,step,'grid');


[~,fi_rand,rad_rand]=cart2sph(x_rand,y_rand,z_rand);
fi_rand=fi_rand*180/pi;
cond = ((fi_rand < 75) & (fi_rand > -75));

r_rand=[x_rand(cond)' y_rand(cond)' z_rand(cond)'];
