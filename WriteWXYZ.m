function WriteWXYZ(lambda,fi,r,w,filename)

in=fopen(filename,'w');

fprintf(in,'%1.13E %1.13E %1.13E %1.13E\n',[lambda(:),fi(:),r(:),w(:)]');

fclose(in);
