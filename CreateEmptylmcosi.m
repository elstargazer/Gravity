function lmcosi=CreateEmptylmcosi(MaxDegree)

Degree=[];
Order=[];

for i=1:MaxDegree+1    
  Degree=[Degree; (i-1)*ones(i,1)];
  Order=[Order; (0:i-1)'];    
end


lmcosi(:,1)=Degree;
lmcosi(:,2)=Order;
lmcosi(:,3)=0;
lmcosi(:,4)=0;