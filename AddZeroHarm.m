function lmcosi=AddZeroHarm(lmcosi,Value)

s=size(lmcosi);
first_line=zeros(1,s(2));
first_line(3)=Value;
lmcosi=[first_line;lmcosi];