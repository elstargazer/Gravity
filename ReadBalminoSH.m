function varargout=ReadBalminoSH(FileName)

if (nargout==1)
% FileName='file_harmo_pot';

in=fopen(FileName,'r');

data=fgets(in);
data=fgets(in);
data=fgets(in);

Rref=str2num(data(1:20));
invf=str2num(data(21:40));
GM=str2num(data(41:60));
omega=str2num(data(61:80));

data=fgets(in);
data=fgets(in);

MaxDegree=str2num(data(17:end));

data=fgets(in);

data=fscanf(in,'%4i%4i   %E %E %E %E %i\n');

fclose(in);

lmcosi=zeros(numel(data)/7,4);

progressbar('Reading','Reordering')

Nu1=numel(data);

for i=1:7:Nu1 
    lmcosi_r(fix(i/7)+1,1)=data(i);
    lmcosi_r(fix(i/7)+1,2)=data(i+1);
    lmcosi_r(fix(i/7)+1,3)=data(i+2);
    lmcosi_r(fix(i/7)+1,4)=data(i+3);
    progressbar(i/Nu1,[]);
end

progressbar(1,[]);

Degree=[];
Order=[];

lmcosi=lmcosi_r;

for i=1:MaxDegree+1    
  Degree=[Degree; (i-1)*ones(i,1)];
  Order=[Order; (0:i-1)'];    
end

Nu2=numel(lmcosi(:,1));

for i=1:Nu2
    
    Index=( (lmcosi_r(:,1)==Degree(i)) & (lmcosi_r(:,2)==Order(i))  );
    lmcosi(i,1)=Degree(i);
    lmcosi(i,2)=Order(i);
    lmcosi(i,3)=lmcosi_r(Index,3);
    lmcosi(i,4)=lmcosi_r(Index,4);  
    
    progressbar(1,i/Nu2);
    
end

progressbar(1,1);

varargout={lmcosi};

elseif (nargout==2)
    
    lmcosi=ReadBalminoSH(FileName);


for i=2:numel(lmcosi(:,1))
    
    C_gt(lmcosi(i,1),lmcosi(i,2)+1)=lmcosi(i,3);
    S_gt(lmcosi(i,1),lmcosi(i,2)+1)=lmcosi(i,4);
    
end

varargout{1}=C_gt;
varargout{2}=S_gt;

end

end

