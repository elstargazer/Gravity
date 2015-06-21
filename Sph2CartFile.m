function Sph2CartFile(InFileName,OutFileName)

data=load(InFileName);

x=data(:,1);%+rand(numel(data(:,1)),1)*50;
y=data(:,2);%+rand(numel(data(:,1)),1)*50;
z=data(:,3);%

[lambda,fi,r]=cart2sph(x,y,z);

lambda=lambda*180/pi;
fi=fi*180/pi;

WriteXYZ(lambda,fi,r,OutFileName);