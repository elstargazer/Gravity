function lmcosi=ReadBalminoSH2(FileName)

in=fopen(FileName,'r');


comment=textscan(in,'%80c\n',1);
data=textscan(in,'%64c\n',1);
data=textscan(in,'%20f64%20f64%20f64%20f64\n',1);
Rref=data{1};
data=textscan(in,'reference date : %f\n',1);
data=textscan(in,'maximal degree : %f\n',1);
MaxDegree=data{1};
data=textscan(in,' l  m dot         cbar                sbar             sigma c      sigma s  lib\n',1);

NumberOfRows=(MaxDegree+1)*(MaxDegree+2)/2;

if (MaxDegree<=999)    
    data=textscan(in,'%3u%3u %21f64%21f64%13f64%13f64 %u\n',NumberOfRows);
elseif (MaxDegree<=9999)
    data=textscan(in,'%4u%4u %21f64%21f64%13f64%13f64 %u\n',NumberOfRows);
else
    data=textscan(in,'%5u%5u %21f64%21f64%13f64%13f64 %u\n',NumberOfRows);
end    

Degree=data{1};
Order=data{2};
C=data{3};
S=data{4};

lmcosi=CreateEmptylmcosi(MaxDegree);

Index=Degree.*(Degree+1)/2+Order+1;

C=C.*Rref;
S=S.*Rref;

lmcosi(Index,3)=C;
lmcosi(Index,4)=S;



fclose(in);

