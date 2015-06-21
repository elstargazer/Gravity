function [dcv dsv]=VerticalDerivatives(m, z, rsquared, Rref,mu,Vmm, Wmm, MaxDeg)

dcv=zeros((MaxDeg-m+1),1);
dsv=zeros((MaxDeg-m+1),1);

    
V0=Vmm;
W0=Wmm;

dcv(1)=V0*mu/Rref;
dsv(1)=W0*mu/Rref;


    
V1=VWVerticalRecursion(m+1,m,z,rsquared,Rref,V0,0);
W1=VWVerticalRecursion(m+1,m,z,rsquared,Rref,W0,0);

dcv(2)=V1*mu/Rref;
dsv(2)=W1*mu/Rref; 

j=2;
    
for n=m+2:MaxDeg        

    Vnew=VWVerticalRecursion(n,m,z,rsquared,Rref,V1,V0);
    Wnew=VWVerticalRecursion(n,m,z,rsquared,Rref,W1,W0);        
        
            
      
    V0=V1;
    W0=W1; 
    V1=Vnew;
    W1=Wnew;
    
    j=j+1;
    
    dcv(j)=V1*mu/Rref;
    dsv(j)=W1*mu/Rref; 
        
end
