function WriteXYZ(lambda,fi,r,filename)

in=fopen(filename,'w');

fprintf(in,'%1.13E %1.13E %1.13E\n',[lambda(:),fi(:),r(:)]');

fclose(in);
