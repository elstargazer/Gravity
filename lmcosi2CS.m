function [C,S]=lmcosi2CS(lmcosi)

MaxDegree=lmcosi(end,1);

lmcosi=AddZeroHarm(lmcosi,1);

for i=1:MaxDegree
    for j=1:i+1
        
        
        C(i,j)=getC(lmcosi,i,j-1);
        S(i,j)=getS(lmcosi,i,j-1);
                
        
    end
end
