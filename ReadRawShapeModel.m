function [lambda,fi,R]=ReadRawShapeModel(FileName,par1,par2,par3,range)

irange=range(1);
jrange=range(2);

in = fopen(FileName, 'r');
R = fread(in,[irange jrange],'uint16','b');
 fclose(in);

R(R==0)=NaN;
 
 
R=par1-par2+(R-1)*par3;


i=1:irange;
j=1:jrange;

i_index=repmat(i',[1 jrange]);
lambda=(i_index-8641)/48;
clear i_index jrange

j_index=repmat(j',[1 irange])';
fi=(4321-j_index)/48;
clear j_index irange

%WriteXYZ(lambda,fi,R,'VestaShapeHamo2.txt');
