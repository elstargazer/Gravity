ccc


Period=[10 9 8 7 6 4.95 4 3];
SurfaceFlattening_n=[0.035459 0.0441583 0.054886 0.0717188 0.0992291 0.1471 0.235441 0.496969];
CoreFlattening_n=[0.0241492 0.0294099 0.0372301 0.0500783 0.0701641  0.099 0.15168 0.263392];


Period_t=2.8:0.02:10.2;


As=3.47647;
Bs=1.31325;

Ac=2.39372;
Bc=1.31325;

SurfaceFlattening_t=(As.*(1-Bs./Period_t./Period_t))./Period_t./Period_t;
CoreFlattening_t=(Ac*(1-Bc./Period_t./Period_t))./Period_t./Period_t;



SurfaceFlattening_nc=interp1(Period,SurfaceFlattening_n,Period_t,'cubic');
CoreFlattening_nc=interp1(Period,CoreFlattening_n,Period_t,'cubic');

%% f plot




figuref = figure('XVisual',''); 
axes1 = axes('Parent',figuref,'FontSize',24);
hold on;

plot(Period_t, SurfaceFlattening_t,'-r','LineWidth',1)
plot(Period_t, CoreFlattening_t,'-b','LineWidth',1)

plot(Period_t, SurfaceFlattening_nc,'--r','LineWidth',1)
plot(Period_t, CoreFlattening_nc,'--b','LineWidth',1)





f2Iv=[0.074,0.090,0.121,0.129,0.131,0.172,0.181,0.178];
f1Iv=[0.052,0.067,0.086,0.090,0.096,0.135,0.134,0.132];
TIv =[6.842,6.000,5.342,5.342,5.000,4.400,4.400,4.400];

% plot(TIv,f1Iv,'o','MarkerSize',7,'MarkerFaceColor','b','MarkerEdgeColor','k');
% plot(TIv,f2Iv,'o','MarkerSize',7,'MarkerFaceColor','r','MarkerEdgeColor','k');
    

xlim([2.8 10.2])
ylim([0 0.5])





xlabel('Rotation period [hr]','FontSize',20);
ylabel('Flattening = (a-c)/a','FontSize',20);

set(gca,'FontSize',20);

set(gcf, 'Units','centimeters', 'Position',[0 0 13 9]*1.5)
set(gcf, 'PaperPositionMode','auto')



data=load('SixOrderHydroList');

plot(data(3,:),data(1,:),'.-r','LineWidth',1);
plot(data(3,:),data(2,:),'.-b','LineWidth',1);

HydrostaticEquilibriumExact2lTest

figure(figuref)


legend({'1st order outer shape','1st order core','Numerical outer shape','Numerical core',...
  '6th order outer shape','6th order core','Exact outer shape','Exact core'})




plot(Period, SurfaceFlattening_n,'or','MarkerSize',5)
plot(Period, CoreFlattening_n,'ob','MarkerSize',5)

CurrentRotationPeriod=5.342;
HydrostaticRotationPeriod=4.95054;
line([CurrentRotationPeriod CurrentRotationPeriod],[0 0.5],'Color','k')
line([HydrostaticRotationPeriod HydrostaticRotationPeriod],[0 0.5],'Color','m')

NorthFlattening=0.1471;
GlobalFlattening=0.1931;

line([2.8 10.2],[NorthFlattening NorthFlattening],'Color','m');
line([2.8 10.2],[GlobalFlattening GlobalFlattening],'Color','k');



box on

print(figuref, '-dpsc', 'HydrostaticEqulibriumFigure.eps');
%% J2 plot

J2current=0.0318;


A=0.564679;


J2_n=[0.00565356, 0.0070497, 0.00878216, 0.0115177, 0.016029, 0.0205784, 0.0233216, 0.0393274, 0.093186];

J2_t= A./(Period_t.^2);

J2_nc=interp1(Period,J2_n,Period_t,'cubic');

figure1 = figure('XVisual',''); 
axes1 = axes('Parent',figure1,'FontSize',24);
hold on;


plot(Period_t, J2_t,'-k','MarkerSize',5);

plot(Period_t,J2_nc,'--k','LineWidth',1)

plot(Period, J2_n,'ok','MarkerSize',5)



legend({'1st order theory','Numerical'})

xlabel('Rotation period [hr]','FontSize',12);
ylabel('J_{2}','FontSize',12);

line([2.8 10.2],[J2current J2current],'Color','k');
line([CurrentRotationPeriod CurrentRotationPeriod],[0 0.5],'Color','k')

set(gca,'FontSize',12);

xlim([2.8 10.2])
ylim([0 0.1])

set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

box on


print(figure1, '-dpsc', 'J2EqulibriumFigure.eps');












