ccc

r1s=261600;
r2s=114200;

filename='/Users/antonermakov/Dawn/grid_search_flattening_075H';

data=load(filename);

a1=data(:,1);
a2=data(:,2);

v1=data(:,4);
v2=data(:,3);

mina1=min(a1(:));
mina2=min(a2(:));

maxa1=max(a1(:));
maxa2=max(a2(:));

a2range=maxa2-mina2;
a1range=maxa1-mina1;

[a1i,a2i]=meshgrid(mina1:a1range/1000:maxa1,mina2:a2range/1000:maxa2);

v1i=griddata(a1,a2,v1,a1i,a2i,'cubic');
v2i=griddata(a1,a2,v2,a1i,a2i,'cubic');

c1i=(r1s.^3)./(a1i.^2);
c2i=(r2s.^3)./(a2i.^2);

f1i=(a1i-c1i)./(a1i);
f2i=(a2i-c2i)./(a2i);



% contour(a1i,a2i,sqrt(v1i.^2+v2i.^2),100)

% contour(a1i,a2i,v2i,100)

% gin=ginput(1)

figure_profile_c=figure();

set(gca,'FontSize',12);
hold on;

v=(sqrt(v1i.^2+v2i.^2));

contour(f1i,f2i,v,30,'Color','k');
pcolor(f1i,f2i,v); 
shading interp;

xlabel('Outer shape flattering','FontSize',12);
ylabel('Core flattening','FontSize',12);


cb=colorbar

set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

ylabel(cb,'Velocity [m/s]','FontSize',12);


%% Find minimum mumerically

Vel = @(X)VelInterp(a1,a2,v1,v2,X(1),X(2));



[xm,fval] = fminsearch(Vel,[mean(a1(:)) mean(a2(:))]);

a1min=xm(1);
a2min=xm(2);


c1min=(r1s.^3)./(a1min.^2);
c2min=(r2s.^3)./(a2min.^2);

f1min=(a1min-c1min)./(a1min);
f2min=(a2min-c2min)./(a2min);


plot(f1min,f2min,'g*','MarkerSize',4);

%% Printing

PrintWhite('FindingHydrostaticFlattenings.eps');














