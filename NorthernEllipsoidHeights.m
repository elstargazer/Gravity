ccc

%% Load shape model

FileName='VestaHASTALAVESTAshape_geoc_elev_3min_mean_g.grd';
%  FileName='VestaHASTALAVESTAshape_geoc_elev_3deg_mean_g.grd';

[ri_shape,fii,lambdai]=ReadRawGrid(FileName);

s=size(ri_shape);

[x,y,z]=sph2cart(lambdai,fii,ri_shape);

x=x(:);
y=y(:);
z=z(:);



x=x_raw(:);
y=y_raw(:);
z=z_raw(:);

% center, radii, evecs

%% Load orientation
data=load('NorthEllipsoidOrientation.mat');


pack_size=10000;

i=[1:pack_size:numel(x) numel(x) ];


altback=zeros(size(x));

progressbar('Progress');

for j=1:numel(i)-1

    
index=i(j):i(j+1);

alt_back(index)=BestFitEllipsoidHeight(x(index),y(index),z(index),center,radii,evecs);
progressbar(j/numel(i))

end

progressbar(1);
% 
% alt_back2=reshape(alt_back,s);
% 
% MapRadialGrid(flipud(ri_shape));




WriteXYZ(lambdai*180/pi,fii*180/pi,alt_back,'NorthernHeight.xyz');


open NorthAltHist.fig

h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;

% NumberOfSlopes=numel(x_raw_cut); 
% y_max=max(max(get(gca,'YTick'))); 
% max_ytick=y_max/NumberOfSlopes*100;
% ylabels = get(gca, 'YTickLabel');
% ylabels = linspace(0,100,length(ylabels));
% set(gca,'YTickLabel',(ylabels*max_ytick/100));
% set(gca,'YTickLabel',fix(10*ylabels*max_ytick/100)/10);


hist(alt_back(fix(numel(x_raw_cut)*rand(1,numel(x_raw_cut)))),100,'r');

% NumberOfSlopes=numel(alt_back); 
% y_max=max(max(get(gca,'YTick'))); 
% max_ytick=y_max/NumberOfSlopes*100;
% ylabels = get(gca, 'YTickLabel');
% ylabels = linspace(0,100,length(ylabels));
% set(gca,'YTickLabel',(ylabels*max_ytick/100));
% set(gca,'YTickLabel',fix(10*ylabels*max_ytick/100)/10);


h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);