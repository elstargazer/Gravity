function [d,CSness,Degree,Order]=GravityPotentialDerivatives(mu,Rref,lmcosi,x,y,z)

% CSness=zeros(437,1);
% Degree=zeros(437,1);
% Order=zeros(437,1);

CSness=[];
Degree=[];
Order=[];

if (lmcosi(1,1)~=0)
s=size(lmcosi);

first_line=zeros(1,s(2));
first_line(3)=1;

lmcosi=[first_line;lmcosi];

end

MaxDeg=lmcosi(end,1);
r=sqrt(x.*x+y.*y+z.*z);
    
Q=0;
rsquared=r.*r;
    
V00=Rref./r;
W00=0;

V0=V00;
W0=W00;

%% Zeroth order (zonal terms)
    
[dcv, ~]=VerticalDerivatives(0,z,rsquared,Rref,mu,V0,W0,MaxDeg);    
    
d=dcv(3:end);

CSness=[CSness; ones(size(dcv(3:end)))];
Degree=[Degree; (2:MaxDeg)'];
Order=[Order; zeros(size(dcv(3:end)))];



%% First order
m=1;

    [V1,W1]=VWDiagonalRecursion(m,x,y,rsquared,Rref,V0,W0);
    
    [dcv, dsv]=VerticalDerivatives(m,z,rsquared,Rref,mu,V1,W1,MaxDeg);
        
    W0=W1;
    V0=V1;  
    
    d=[d; dcv(2:end); dsv(2:end)];
    
    CSness=[CSness; ones(size(dcv(2:end))); zeros(size(dsv(2:end)))];
    Degree=[Degree; (2:MaxDeg)'; (2:MaxDeg)'];
    Order=[Order; ones(size(dcv(2:end))); ones(size(dsv(2:end)))];
 
%% Second and higher orders

for m=2:MaxDeg-2   

    V1=VWDiagonalRecursion(m,x,y,rsquared,Rref,V0,W0);
    W1=VWDiagonalRecursion(m,x,y,rsquared,Rref,W0,-V0);
    
    [dcv, dsv]=VerticalDerivatives(m,z,rsquared,Rref,mu,V1,W1,MaxDeg);
        
    W0=W1;
    V0=V1;  
    
    d=[d; dcv; dsv];
    
    CSness=[CSness; ones(size(dcv)); zeros(size(dsv))];
    Degree=[Degree; (m:MaxDeg)'; (m:MaxDeg)'];
    Order=[Order; m*ones(size(dcv)); m*ones(size(dsv))];
    
end

%% Penultimate order
    
[V1,W1]=VWDiagonalRecursion(MaxDeg-1,x,y,rsquared,Rref,V0,W0);
    
V0=V1;
W0=W1;

d=[d; V1*mu/Rref; W1*mu/Rref];

CSness=[CSness;1;0 ];
Degree=[Degree; MaxDeg-1; MaxDeg-1];
Order=[Order; MaxDeg-1; MaxDeg-1];


V1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,Rref,W0,0); 

d=[d; V1*mu/Rref; W1*mu/Rref];

CSness=[CSness;1;0 ];
Degree=[Degree; MaxDeg; MaxDeg];
Order=[Order; MaxDeg-1; MaxDeg-1];

%% Last order
    
Vlast=VWDiagonalRecursion(MaxDeg,x,y,rsquared,Rref,V0,W0);
Wlast=VWDiagonalRecursion(MaxDeg,x,y,rsquared,Rref,W0,-V0);

d=[d; Vlast*mu/Rref; Wlast*mu/Rref];

CSness=[CSness;1;0];
Degree=[Degree; MaxDeg; MaxDeg];
Order=[Order; MaxDeg; MaxDeg];

d=d./NormCoef(Degree,Order);





