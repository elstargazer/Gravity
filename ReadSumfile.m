function [x,y,z,t]=ReadSumfile(filename)

in_file=fopen(filename,'r');

str=fgetl(in_file);

time_string=fgetl(in_file);
t=[time_string ' UTC'];

str=fgetl(in_file);
str=fgetl(in_file);
str=fgetl(in_file);

C=textscan(str,'%f %f %f%*[^\n]');

x=C{1};
y=C{2};
z=C{3};


fclose(in_file);


    