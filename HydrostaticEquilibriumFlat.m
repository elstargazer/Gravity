ccc
%% Input parameters


GM=17.2883080293e9;
G=6.67e-11;
M=GM/G;

V=0.74886e8*1e9;

Rref=265000;

h=40000; % crustal thickness

J2obs=0.03177947410998;
J2outer=0.0351861;


router=(V/(4/3*pi))^(1/3);

T=4.88:0.02:4.98;
% W=2*pi/(T*3600);

rcore=85000:1000:220000;

fouter_obs=0.147;

fcore0=0.1;
fouter0=0.1;

% rhocore=4000:25:8000;


%% Plotting settings
fig1=figure; hold on; 
ax1=gca;
set(ax1,'FontSize',12);

data=load('ListJ2.txt');


xlim([85 220]);
ylim([3000 8000]);

xlabel('Core size [km]','FontSize',12);
ylabel('Core density [kg/m^3]','FontSize',12);

rhocore0=5000;


%% Grid of core radii and densities
rhocoreg=3000:100:8000;

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);

rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));


% progressbar(0)

for i=1:numel(rhocorei)

    [fhi(i,:),fval]=HydrostaticStateExact2l(router,rcorei(i),4.9,rhoouteri(i),rhocorei(i),0.1,0.1);
    fouter0=fhi(i,1);
    fcore0=fhi(i,2);
    
    progressbar(i/numel(rhocorei));
    

end

progressbar(1);


fcorei=fhi(:,2);

fcorei=reshape(fcorei,size(rhocorei));


% surf(rcorei/1000,rhocorei,fcorei); shading interp
% caxis([0.08 0.16]);

% alpha(0.1);
%% Contour mantle density
levels=[0:500:2000 2250:250:5000];

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri,levels,'Color',[0.4 0.4 0.4]);
text_handle = clabel(C,h,'manual');
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7]);



%% Contour core flattening

levels=0:0.01:0.2
[C,h]=contour(rcorei/1000,rhocorei,fcorei,levels,'Color',[0 1 0]);
text_handle = clabel(C,h,'manual');

% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7]);

%% Plotting J2 solution

% rcore_m=85:1:230;
% 
% plot(rcore_m,data(:,1),'r')
% plot(rcore_m,data(:,2),'r')


% f1obs=0.19;
f2max=0.163;
f2min=0.098;

rhocore_j2soln_1=CoreDensityFromJ2Mass(M,rcore,f2max,Rref,J2outer,J2obs,V);
rhocore_j2soln_2=CoreDensityFromJ2Mass(M,rcore,f2min,Rref,J2outer,J2obs,V);

% 


plot(ax1,rcore/1000,rhocore_j2soln_1,'r')
plot(ax1,rcore/1000,rhocore_j2soln_2,'r')

X=[rcore/1000,fliplr(rcore/1000)];                %#create continuous x value array for plotting
Y=[rhocore_j2soln_1 fliplr(rhocore_j2soln_2)];              %#create y values for out and then back
fill(X,Y,'r'); 

%% Main loop

Tol=0.0001;
Eps=Tol+1;

progressbar(0, 0)

for k=1:numel(T)
    

fcore0=0.1;
fouter0=0.1;


    for i=1:numel(rcore)
        
   
    
    rhocore_h(i)=fzero(@(rhocore) FunFlat(rhocore,M,rcore(i),router,T(k),fouter_obs),rhocore0);
    
    rhoouter=-(3*M-4*pi*((rcore(i)).^3).*rhocore_h(i))./(4*pi*(rcore(i).^3)-4*pi*(router^3));
    
    [fh,d2]=HydrostaticStateExact2l(router,rcore(i),T(k),rhoouter,rhocore_h(i),fouter0,fcore0);
    
    rhocore_h2(i)=CoreDensityFromJ2Mass(M,rcore(i),fh(2),Rref,J2outer,J2obs,V);
    
    
        while (Eps>Tol)
            
        rhoouter=-(3*M-4*pi*((rcore(i)).^3).*rhocore_h2(i))./(4*pi*(rcore(i).^3)-4*pi*(router^3));
            
        fhold=fh;
        
        [fh,d2]=HydrostaticStateExact2l(router,rcore(i),T(k),rhoouter,rhocore_h2(i),fh(1),fh(2));
        
        rhocore_h2(i)=CoreDensityFromJ2Mass(M,rcore(i),fh(2),Rref,J2outer,J2obs,V);
        
        Eps=max(abs(fhold-fh));
        
        end
     
    rhocore0=rhocore_h(i);
    progressbar(k/numel(T), i/numel(rcore));
    
    fouter0=fh(1);
    fcore0=fh(2);
    
    
    end
    
    

    
    % plot solution for given T
    
    plot(ax1,rcore/1000,rhocore_h,'-k','LineWidth',1);
    plot(ax1,rcore/1000,rhocore_h2,'-b','LineWidth',1);
    
    progressbar(k/numel(T), i/numel(rcore));
    
    rhocore0=rhocore_h(1);


end

progressbar(1,1);

%% Plotting




set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

box on


print(fig1, '-dpsc', 'HydrostaticEquilibriumExactSolutionJ2flat.eps');















