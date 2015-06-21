clear all;
close all;

ShapeFileNameStereo='stereo_shape.xyz';
GridFileNameStereo='stereo_grid.xyz';


data1=load(ShapeFileNameStereo);
data2=load(GridFileNameStereo);

size(data1)
size(data2)

xs=data1(:,1);
ys=data1(:,2);
zs=data1(:,3);

xg=data2(:,1);
yg=data2(:,2);
zg=data2(:,3);

condition=zg==zg(1);
zg(condition)=NaN;
condition=isnan(zg);
%%
minxg=min(min(xg));
maxxg=max(max(xg));
minyg=min(min(yg));
maxyg=max(max(yg));

xg_range=maxxg-minxg;
yg_range=maxyg-minyg;

xg_range
yg_range


%%
minxs=min(min(xs));
maxxs=max(max(xs));
minys=min(min(ys));
maxys=max(max(ys));

xs_range=maxxs-minxs;
ys_range=maxys-minys;

xs_range
ys_range

%%
[minxg minyg]
[minxs minys]

[maxyg maxyg]
[maxys maxys]


zgi=griddata(xs,ys,zs,xg,yg,'nearest');
zg(condition)=0;
zgi(~condition)=0;
zg=zg+zgi;

sum(isnan(zg))/numel(zg)
