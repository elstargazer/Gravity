function [n,Z]=SphericalHarmonicAdmittance(lmcosi_grav,lmcosi_shape,mu,Rref)

R0 = lmcosi_shape(1,3);

MaxDegree=min(lmcosi_grav(end,1),lmcosi_shape(end,1));

deg = lmcosi_grav(:,1); 
lmcosi_grav(:,3)=lmcosi_grav(:,3).*(deg+1)*mu/(Rref^2).*1e5;
lmcosi_grav(:,4)=lmcosi_grav(:,4).*(deg+1)*mu/(Rref^2).*1e5;

n = 2:MaxDegree;
Z=zeros(1,numel(n));

for i=1:numel(n)   
    i1=(n(i))*(n(i)+1)/2+1;
    i2=(n(i)+1)*(n(i)+2)/2;   
  
    Z(i)=sum(lmcosi_grav(i1:i2,3).*lmcosi_shape(i1:i2,3)+....
             lmcosi_grav(i1:i2,4).*lmcosi_shape(i1:i2,4))./...
         sum(lmcosi_shape(i1:i2,3).^2+lmcosi_shape(i1:i2,4).^2);       
end

Z=Z*1000;

%((R0/Rref).^deg)