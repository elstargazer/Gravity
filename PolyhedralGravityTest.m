ccc

G=6.67384e-11;

%% Create shape
NTess=0;
TR=IcosahedronMesh;
TR_2=SubdivideSphericalMesh(TR,NTess); 

FV=TR_2.Triangulation;
x1=TR_2.X(:,1)';
y1=TR_2.X(:,2)';
z1=TR_2.X(:,3)';

[lambdai,fii,ri]=cart2sph(x1,y1,z1);

s=size(FV);
Nelems=s(1);

rho=1*ones(1,Nelems);

n=1;
rho_m=rho*0;
rho_m(n)=rho_m(n)+1;

%% Create reference surface
Rref=1.01;
FiStep=2;
LambdaStep=2;
[lambdare,fire]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );

xre=Rref*cos(lambdare).*cos(fire);
yre=Rref*sin(lambdare).*cos(fire);
zre=Rref*sin(fire);

tic
U1  =GravityPotentialTri93   (x1,y1,z1,FV,xre,yre,zre,rho(1));
toc

% tic
% U1_m=GravityPotentialTriDen93(x1,y1,z1,FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho_m );
% toc

U1_m=0;

progressbar(0);

for i=1:Nelems
    
    rho_m=rho*0;
    rho_m(i)=1;
    U1_m=U1_m+GravityPotentialTriDen93(x1,y1,z1,FV,xre,yre,zre,rho_m );
    progressbar(i/Nelems);
    
end

progressbar(1);
   

% U1der=GravityPotentialTriDenDer93(x1,y1,z1,FV(:,1)',FV(:,2)',FV(:,3)',xre(:),yre(:),zre(:));
% U1der_part=reshape(U1der(:,n),size(xre));


U_diff=U1_m-U1;

% U2=GravityPotentialTriDen(x1,y1,z1,FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho);
% U2_m=GravityPotentialTriDen(x1,y1,z1,FV(:,1)',FV(:,2)',FV(:,3)',xre,yre,zre,rho_m);


AGUaxes; hold on;
pcolorm(fire,lambdare,U_diff); shading interp;
colorbar('FontSize',20);
plotm(fii,lambdai,'wo','MarkerFaceColor','w','MarkerSize',5);

plotm(fii(FV(n,1)),lambdai(FV(n,1)),'ko','MarkerFaceColor','k','MarkerSize',8);
plotm(fii(FV(n,2)),lambdai(FV(n,2)),'ko','MarkerFaceColor','k','MarkerSize',8);
plotm(fii(FV(n,3)),lambdai(FV(n,3)),'ko','MarkerFaceColor','k','MarkerSize',8);


% AGUaxes
% pcolorm(fire,lambdare,U2); shading interp;
% colorbar('FontSize',20);
% 
% 
% AGUaxes
% pcolorm(fire,lambdare,U2-U1); shading interp;
% colorbar('FontSize',20);


% quiver3(xr,yr,zr,-ax,-ay,-az,'r')

% Vsph=4/3*pi
% V=TetrahedronBobyVolume(x1,y1,z1,FV(:,1)',FV(:,2)',FV(:,3)');
% (Vsph-V)/Vsph
% 
% Upoint=G*V*rho(1)/Rref
% U1mean=mean(U1(:))
% U1mean_m=mean(U1_m(:))

% U2mean=mean(U2(:));