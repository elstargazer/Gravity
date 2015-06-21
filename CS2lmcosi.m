function lmcosi=CS2lmcosi(C,S)

s=size(C);

MaxDegree=s(1);

for i=1:MaxDegree
    
    i_start=(i)*(i+1)/2+1;
    i_end=(i+1)*(i+2)/2;
    
    
    lmcosi(i_start:i_end,1)=i;
    lmcosi(i_start:i_end,2)=0:i;
    lmcosi(i_start:i_end,3)=C(i,1:i+1);
    lmcosi(i_start:i_end,4)=S(i,1:i+1);    
    
end

lmcosi=lmcosi(2:end,:);