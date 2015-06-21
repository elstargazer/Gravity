function U=GravityPotential2(mu,Rref,lmcosi,x,y,z)


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
    
Q=Q+VerticalSum(0,z,rsquared,Rref,V0,W0,MaxDeg,lmcosi);    
    
    
for m=1:MaxDeg-2   

    [V1,W1]=VWDiagonalRecursion(m,x,y,rsquared,Rref,V0,W0);
    
    Q=Q+VerticalSum(m,z,rsquared,Rref,V1,W1,MaxDeg,lmcosi);
        
    W0=W1;
    V0=V1;        
end
    
[V1,W1]=VWDiagonalRecursion(MaxDeg-1,x,y,rsquared,Rref,V0,W0);
    
V0=V1;
W0=W1;
    
[Vlast,Wlast]=VWDiagonalRecursion(MaxDeg,x,y,rsquared,Rref,V0,W0);
    
V1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,Rref,W0,0);  
    
Q=Q+(getC(lmcosi,MaxDeg-1,MaxDeg-1).*V0+getS(lmcosi,MaxDeg-1,MaxDeg-1).*W0)./NormCoef(MaxDeg-1,MaxDeg-1);
Q=Q+(getC(lmcosi,MaxDeg,MaxDeg-1).*V1+getS(lmcosi,MaxDeg,MaxDeg-1).*W1)./NormCoef(MaxDeg,MaxDeg-1);

Q=Q+(getC(lmcosi,MaxDeg,MaxDeg).*Vlast+getS(lmcosi,MaxDeg,MaxDeg).*Wlast)./NormCoef(MaxDeg,MaxDeg); 
 
% a=(getC(lmcosi,MaxDeg,MaxDeg)*Vlast+getS(lmcosi,MaxDeg,MaxDeg)*Wlast)/NormCoef(MaxDeg,MaxDeg);
% 
% disp('a=');
% a*mu/Rref


U=mu.*Q./Rref;
