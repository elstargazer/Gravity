function S=CrossPower(lmcosi1,lmcosi2)

MaxDegree=min(lmcosi1(end,1),lmcosi2(end,1));

S=zeros(1,MaxDegree);



for i=3:MaxDegree+1;
    
    i_start=(i-1)*(i)/2+1;
    i_end=(i)*(i+1)/2;   
    
    
    
    S(i-1)=sum(lmcosi1(i_start:i_end,3).*lmcosi2(i_start:i_end,3)+lmcosi1(i_start:i_end,4).*lmcosi2(i_start:i_end,4));
         
end