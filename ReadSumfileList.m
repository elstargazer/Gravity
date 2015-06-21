function [x_out,y_out,z_out,t_out]=ReadSumfileList(folder,filelist)

str=' ';

in_list=fopen(filelist,'r');


str=fgetl(in_list);

x_out=[];
y_out=[];
z_out=[];
t_out=[];


while (str~=-1)
    
    [x,y,z,t]=ReadSumfile([folder '/' str]);
    x_out=[x_out; x];
    y_out=[y_out; y];
    z_out=[z_out; z];
    t_out=[t_out; t];
    
    str=fgetl(in_list);

end
    
fclose(in_list);
    
    