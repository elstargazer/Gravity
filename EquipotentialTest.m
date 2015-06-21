ccc


NTess=0;
 
TR=IcosahedronMesh;
 
TR_2=SubdivideSphericalMesh(TR,NTess); 
 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), 
% axis equal; 
% hold on; 
% axis([-5 5 -5 5 -5 5]/2);
 
 
FV=TR_2.Triangulation;

T=20;
w=2*pi/(T*3600);

f=0.032:0.0001:0.049;

R=100000;

for i=1:numel(f)

a=R/((1 - f(i))^(1/3));
c=(R - f(i)*R)/((1 - f(i))^(1/3));
  
 
 
xt=a*TR_2.X(:,1)';
yt=a*TR_2.X(:,2)';
zt=c*TR_2.X(:,3)';

 
[U,~,~,~]=GravityPotentialAccelerationTri(xt,yt,zt,FV(:,1)',FV(:,2)',FV(:,3)',xt*1.00001,yt*1.00001,zt*1.00001,1000);
% [U2,~,~,~]=GravityPotentialAccelerationTri(xt,yt,zt,FV(:,1)',FV(:,2)',FV(:,3)',xt*1.00001,yt*1.00001,zt*1.00001,1000);

 
T=+0.5*(w*w).*(xt.*xt+yt.*yt);
V=U+T;



U_std(i)=std(V(:));


end


figure; hold on;
set(gca,'FontSize',20);

plot(f,U_std);
