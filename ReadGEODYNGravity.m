function [lmcosi,mu,Rref]=ReadGEODYNGravity(FileName)

in=fopen(FileName,'r');

str=fgets(in);

str=fgets(in);

MaxDegree=str2num(str(16:17));
mu=str2num(str(25:44));
Rref=str2num(str(45:59));


lmcosi=CreateEmptylmcosi(MaxDegree);


lmcosi(1,3)=1;
lmcosi(2,3)=0;
lmcosi(3,3)=0;
lmcosi(3,4)=0;

str=fgets(in);

while (str~=-1)
    
    
n=str2num(str(16:17));
m=str2num(str(19:20));


row=(lmcosi(:,1)==n) & (lmcosi(:,2)==m);

if (str(1:7)=='GCOEFC1')
    lmcosi(row,3)=str2num(str(25:44));
    lmcosi(row,5)=str2num(str(60:72));
elseif  (str(1:7)=='GCOEFS1')
    lmcosi(row,4)=str2num(str(25:44));
    lmcosi(row,6)=str2num(str(60:72));
end
    
    
    
str=fgets(in);

end





fclose(in)