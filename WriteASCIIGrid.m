function WriteASCIIGrid(FileName,r)


% FiStep=.10;
% LambdaStep=.10;

s=size(r);

N1=s(1);
N2=s(2);

% FiStep=1/20;
% LambdaStep=1/20;
% 
% 
% N1=180/FiStep;
% N2=360/LambdaStep;


in=fopen(FileName,'w');

format='   ';

format_number='%16.9f';

progressbar(0)

for i=1:N2-1
    format=[format format_number '  '];
    
    progressbar(i/(N2-1));
end

progressbar(1);

format=[format format_number '\n'];

progressbar(0)

tic

for i=1:(N1)
    
    fprintf(in,format,r(i,:));
    
    progressbar(i/N1)
    
end

time_el=toc;

time_el


progressbar(1)

fclose(in);