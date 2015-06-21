function q=IsotropicRatio(lmcosi1,lmcosi2)

MaxDegree=min(lmcosi1(end,1),lmcosi2(end,1));

q=zeros(1,MaxDegree);

n=lmcosi1(:,1);
m=lmcosi2(:,2);
c1=lmcosi1(:,3);
s1=lmcosi1(:,4);
c2=lmcosi2(:,3);
s2=lmcosi2(:,4);

for i=2:MaxDegree+1;
    
    i1=(i-1)*(i)/2;
    i2=(i)*(i+1)/2-1;   
    
    a=i1:i2;    
    
    mag=c1(a).*c2(a)+s1(a).*s2(a);
    
    st=sum((n(a).*(n(a)+1)-(2.*n(a)+1).*m(a)/2).*mag);
    sb=sum((2*n(a)+1).*m(a)/2.*mag);
    
    q(i-1)=st/sb;  
end