function [ax,ay,az]=GravityAcceleration(mu,Rref,lmcosi,x,y,z)

MaxDegree=lmcosi(end,1);

if (lmcosi(1,1)~=0)
s=size(lmcosi);

first_line=zeros(1,s(2));
first_line(3)=1;

lmcosi=[first_line;lmcosi];

end
s=size(lmcosi);
MaxDeg=lmcosi(end,1);
lmcosi=[lmcosi; zeros(MaxDegree+2,s(2))];

s=size(x);

ax=zeros(s);
ay=zeros(s);
az=zeros(s);

r=sqrt(x.*x+y.*y+z.*z);
    
factor=mu/Rref/Rref;
rsquared=r.*r;
    
V00=Rref./r;
W00=0;
    
V0=V00;
W0=W00;

%% Zonal terms

m=0;
n=1;

V1=VWVerticalRecursion(n,m,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(n,m,z,rsquared,Rref,W0,0);  
%
az=az+((n-m)*(-getC(lmcosi,n-1,m).*V1-getS(lmcosi,n-1,m).*W1));
%

for n=m+2:MaxDeg+1
    
V=VWVerticalRecursion(n,m,z,rsquared,Rref,V1,V0);
W=VWVerticalRecursion(n,m,z,rsquared,Rref,W1,W0);  

%    
az=az+((n-m)*(-getC(lmcosi,n-1,m).*V-getS(lmcosi,n-1,m).*W))/NormCoef(n-1,m);

ax=ax+0.5*(factorial(n-m)/factorial(n-m-2)*(getC(lmcosi,n-1,m+1)*V+getS(lmcosi,n-1,m+1)*W))/NormCoef(n-1,m+1);
ay=ay+0.5*(factorial(n-m)/factorial(n-m-2)*(-getC(lmcosi,n-1,m+1)*W+getS(lmcosi,n-1,m+1)*V))/NormCoef(n-1,m+1);

%
V0=V1;
W0=W1;

V1=V;
W1=W;

end

%% First order terms

m=1;
n=m;
[Vsect,Wsect]=VWDiagonalRecursion(m,x,y,rsquared,Rref,V00,W00);

V0=Vsect;
W0=Wsect;

%
ax=ax+(-getC(lmcosi,0,m-1)*V0);
ay=ay+(-getC(lmcosi,0,m-1)*W0);
%

n=m+1;

V1=VWVerticalRecursion(n,m,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(n,m,z,rsquared,Rref,W0,0); 

% 
az=az+((n-m)*(-getC(lmcosi,n-1,m).*V1-getS(lmcosi,n-1,m).*W1))/NormCoef(n-1,m);
        
ax=ax+(-getC(lmcosi,n-1,m-1)*V1)/NormCoef(n-1,m-1);
ay=ay+(-getC(lmcosi,n-1,m-1)*W1)/NormCoef(n-1,m-1);
%

for n=m+2:MaxDeg+1
    
       V=VWVerticalRecursion(n,m,z,rsquared,Rref,V1,V0);
       W=VWVerticalRecursion(n,m,z,rsquared,Rref,W1,W0);  
       
       %
       az=az+((n-m)*(-getC(lmcosi,n-1,m).*V-getS(lmcosi,n-1,m).*W))/NormCoef(n-1,m);
        
       ax=ax+(-getC(lmcosi,n-1,m-1)*V)/NormCoef(n-1,m-1);
       ay=ay+(-getC(lmcosi,n-1,m-1)*W)/NormCoef(n-1,m-1);
        
       ax=ax+0.5*(factorial(n-m)/factorial(n-m-2)*(getC(lmcosi,n-1,m+1)*V+getS(lmcosi,n-1,m+1)*W))/NormCoef(n-1,m+1);
       ay=ay+0.5*(factorial(n-m)/factorial(n-m-2)*(-getC(lmcosi,n-1,m+1)*W+getS(lmcosi,n-1,m+1)*V))/NormCoef(n-1,m+1);
       %
       
       V0=V1;
       W0=W1;

       V1=V;
       W1=W;    
end

%% Sectorial and tesseral terms (m>2)

m=2;

[Vsect,Wsect]=VWDiagonalRecursion(m,x,y,rsquared,Rref,Vsect,Wsect);

for m=2:MaxDeg+1

n=m;
V0=Vsect;
W0=Wsect;

%
 ax=ax+0.5*(-getC(lmcosi,m-1,m-1)*V0-getS(lmcosi,m-1,m-1)*W0)/NormCoef(m-1,m-1);
 ay=ay+0.5*(-getC(lmcosi,m-1,m-1)*W0+getS(lmcosi,m-1,m-1)*V0)/NormCoef(m-1,m-1);
% 

n=m+1;

V1=VWVerticalRecursion(n,m,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(n,m,z,rsquared,Rref,W0,0); 

% 
az=az+((n-m)*(-getC(lmcosi,n-1,m).*V1-getS(lmcosi,n-1,m).*W1))/NormCoef(n-1,m);
       
ax=ax+0.5*(-getC(lmcosi,n-1,m-1)*V1-getS(lmcosi,n-1,m-1)*W1)/NormCoef(n-1,m-1);
ay=ay+0.5*(-getC(lmcosi,n-1,m-1)*W1+getS(lmcosi,n-1,m-1)*V1)/NormCoef(n-1,m-1);

%
    for n=m+2:MaxDeg+1
        
       V=VWVerticalRecursion(n,m,z,rsquared,Rref,V1,V0);
       W=VWVerticalRecursion(n,m,z,rsquared,Rref,W1,W0); 
       
       %
       az=az+((n-m)*(-getC(lmcosi,n-1,m).*V-getS(lmcosi,n-1,m).*W))/NormCoef(n-1,m);
       
       ax=ax+0.5*(-getC(lmcosi,n-1,m-1)*V-getS(lmcosi,n-1,m-1)*W)/NormCoef(n-1,m-1);
       ay=ay+0.5*(-getC(lmcosi,n-1,m-1)*W+getS(lmcosi,n-1,m-1)*V)/NormCoef(n-1,m-1);
       
       ax=ax+0.5*(factorial(n-m)/factorial(n-m-2)*(getC(lmcosi,n-1,m+1)*V+getS(lmcosi,n-1,m+1)*W))/NormCoef(n-1,m+1);
       ay=ay+0.5*(factorial(n-m)/factorial(n-m-2)*(-getC(lmcosi,n-1,m+1)*W+getS(lmcosi,n-1,m+1)*V))/NormCoef(n-1,m+1);
       %
       
       V0=V1;
       W0=W1;

       V1=V;
       W1=W;   
       
    end   
    
[Vsect,Wsect]=VWDiagonalRecursion(m+1,x,y,rsquared,Rref,Vsect,Wsect);
       
end

ax=-ax*factor;
ay=-ay*factor;
az=-az*factor;
