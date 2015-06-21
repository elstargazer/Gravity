close all;
clc;
clear all;

filename='ver512q.tab.txt';

in=fopen(filename);

str=fgetl(in);

data=sscanf(str,'%d %d');
N1=data(1);
N2=data(2);

data1=fscanf(in,'%d %f %f %f\n',[N1 4]);

data1=reshape(data1,[4 N1]);
x=data1(2,:);
y=data1(3,:);
z=data1(4,:);

data2=fscanf(in,'%d %d %d %d\n',[N2 4]);

data2=reshape(data2,[4 N2]);
Ver1=data2(2,:);
Ver2=data2(3,:);
Ver3=data2(4,:);


fclose(in);

figure('color','k','Position',[1 1 2560 1440]);

plot3(max(x(:)),max(y(:)),max(z(:)),'k.','MarkerSize',1);
plot3(min(x(:)),min(y(:)),min(z(:)),'k.','MarkerSize',1);

h=trisurf([Ver1; Ver2; Ver3]',x,y,z,z*0); shading interp

lighting phong
shading interp
axis equal tight off
box off
hold on;
light_handle=light('Style','infinite');
set(h,'BackFaceLighting','lit');
colormap gray;
material([0 0.9 0]);

% clear data1 data2
% clear 