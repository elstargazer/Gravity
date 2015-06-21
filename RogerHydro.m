ccc

data=load('grid_search_Sch4.txt');

ao=data(:,1)/1000;
ai=data(:,2)/1000;
v1=data(:,3);
v2=data(:,4);

v=sqrt(v1.*v1+v2.*v2);

aomin=min(ao);
aomax=max(ao);
aimax=max(ai);
aimin=min(ai);

step=0.005;

[aoi,aii]=meshgrid(aomin:step:aomax,aimin:step:aimax);

xs=273263;
ys=216524;


%% 



%% Cubic spline interpolation
vi=griddata(ao,ai,v,aoi,aii,'cubic');

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

% x=(273:0.15:274)-0.58;
% y=215:0.15:219;
% 
% mx=mean(x);
% my=mean(y);
% 
% [xi,yi]=meshgrid(x-mx,y-my);
% 
% xi2=xi(:);
% yi2=yi(:);
% 
% theta=43/180*pi;
% 
% 
% xi3=xi2*cos(theta)-yi2*sin(theta);
% yi3=xi2*sin(theta)+yi2*cos(theta);
% 
% xi3=reshape(xi3+mx,size(xi));
% yi3=reshape(yi3+my,size(yi));
% plot(xi3,yi3,'*k')
% 
% 
% WriteXYZ(xi3(:),yi3(:),yi3(:),'grid_search_rot.txt');
% 
%  data=load('grid_search_rot.txt');
%  plot(data(:,1),data(:,2),'.');
% for i=1:numel(aoi)
%     
%     text(ao(i),ai(i),num2str(v(i)*1e16,'%4.2u'),'FontSize',15);
%       
%     
%     
% end

index=find(v==min(v));

plot(ao(index),ai(index),'*b','MarkerSize',20);


WriteXYZ(ao(:),ai(:),v(:),'grid_search_check.txt');

xs-aoi(index)*1000
ys-aii(index)*1000


%% Polynomial fit

coeffs = fit2dPolySVD( ao, ai, (v),2 );


aoif=aoi(:);
aiif=aii(:);


vbar=eval2dPoly(aoif,aiif,coeffs);


aoif=reshape(aoif,size(aoi));
aiif=reshape(aiif,size(aii));
vbar=reshape(vbar,size(aoi));

index=find(vbar==min(vbar(:)));

figure; hold on;

% plot(ao,ai,'ok')
contour(aoif,aiif,vbar,400);

set(gca,'FontSize',20);
xlabel('Outer shape equatorial radius [km]','FontSize',20);
ylabel('Core shape equatorial radius [km]','FontSize',20);

axis equal

xlim([aomin aomax])
ylim([aimin aimax])

plot(aoif(index),aiif(index),'*k');
plot(xs/1000,ys/1000,'*g');


%%
figure; hold on;

vbar2=eval2dPoly(ao,ai,coeffs);

dvi=griddata(ao,ai,(v-vbar2)./v,aoi,aii,'cubic');

contour(aoi,aii,dvi*100,50);

colorbar

%%

r1=210.570;
r2=265.301;

f1S=0.08024799538136373;
f2S=0.08488811613005487;


f1if=(aiif-(r1.^3)./(aiif.^2))./(aiif);
f2if=(aoif-(r2.^3)./(aoif.^2))./(aoif);


figure; hold on;
set(gca,'FontSize',20);


contour(f1if,f2if,log10(vbar),100)

xlabel('Core flattening','FontSize',20);
ylabel('Outer flattening','FontSize',20);



plot(f1if(index),f2if(index),'k*');

plot(f1S,f2S,'g*');


axis equal



xlim([min(f1if(:)) max(f1if(:))])
ylim([min(f2if(:)) max(f2if(:))])

[f1S f1if(index)]
[f2S f2if(index)]


df1=(f1S-f1if(index))/f1S*100
df2=(f2S-f2if(index))/f2S*100






















