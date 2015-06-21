function Qv=VerticalSum(m, z, rsquared, Rref,Vmm, Wmm, MaxDeg,lmcosi)

Qv=0;
    
V0=Vmm;
W0=Wmm;
    
V1=VWVerticalRecursion(m+1,m,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(m+1,m,z,rsquared,Rref,W0,0);       
    
Qv=Qv+(getC(lmcosi,m,m).*V0+getS(lmcosi,m,m).*W0)./NormCoef(m,m);    
Qv=Qv+(getC(lmcosi,m+1,m).*V1+getS(lmcosi,m+1,m).*W1)./NormCoef(m+1,m);  
    
for n=m+2:MaxDeg        

    Vnew=VWVerticalRecursion(n,m,z,rsquared,Rref,V1,V0);
    Wnew=VWVerticalRecursion(n,m,z,rsquared,Rref,W1,W0);        
        
    Qv=Qv+(getC(lmcosi,n,m).*Vnew+getS(lmcosi,n,m).*Wnew)./NormCoef(n,m);
        
      
    V0=V1;
    W0=W1; 
    V1=Vnew;
    W1=Wnew;
        
end