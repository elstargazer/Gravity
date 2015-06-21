ccc

data=load('grid_search_Sch4.txt');

ao=data(:,1)/1000;
ai=data(:,2)/1000;
v1=data(:,3);
v2=data(:,3);

v=sqrt(v1.*v1+v2.*v2);

aomin=min(ao);
aomax=max(ao);
aimax=max(ai);
aimin=min(ai);

step=0.005;

[aoi,aii]=meshgrid(aomin:step:aomax,aimin:step:aimax);

xs=273263;
ys=216524;


%% Cubic spline interpolation
F = TriScatteredInterp(ao,ai,v);


vi = F(aoi,aii);

figure; hold on;

plot(ao,ai,'ok')
contour(aoi,aii,vi,100);
%  surf(aoi,aii,vi)

% surf(aoi,aii,v)

set(gca,'FontSize',20);
xlabel('Outer shape equatorial radius [km]','FontSize',20);
ylabel('Core shape equatorial radius [km]','FontSize',20);


% axis equal

xlim([aomin aomax])
ylim([aimin aimax])


index=find(vi==min(vi(:)));

plot(aoi(index),aii(index),'*g','MarkerSize',10);

plot(xs/1000,ys/1000,'k*');

index=find(v==min(v));

plot(ao(index),ai(index),'*b','MarkerSize',20);