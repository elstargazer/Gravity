function [lmcosi,Rref] = ReadRyanSH(filename)

in = fopen(filename);

for i=1:5
    fgetl(in);
end

C = textscan(in,'  OBR   = %f\',1);
Rref = C{1}*1000;

fgetl(in);

data  = textscan(in,'OBAJ(%d8) =%f,');
nz = double(data{1});
cz = -data{2};

data  = textscan(in,'OBAC(%d8,%d8) =%f,  OBAS(%d8,%d8) =%f,\n');
n  = double(data{1});
m  = double(data{2});
c  = data{3};
s  = data{6};

nmax = max(nz);

lmcosi = CreateEmptylmcosi(nmax);

indz = (nz+2).*(nz+1)/2-nz;
lmcosi(indz,3) = cz;

ind  = (n+2).*(n+1)/2 - n + m;

lmcosi(ind,3) = c;
lmcosi(ind,4) = s;

fclose(in);

