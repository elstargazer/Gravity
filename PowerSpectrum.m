function [k,sdl]=PowerSpectrum(lmcosi)
% as in Turcotte 1987

s=size(lmcosi);

MaxDegree=lmcosi(s(1),1);

Rref=lmcosi(1,3)/1000;

lmcosi(:,3)=lmcosi(:,3)/Rref;
lmcosi(:,4)=lmcosi(:,4)/Rref;


for i=1:MaxDegree
    
   a1=(i)*(i+1)/2+1;
   a2=(i+1)*(i+2)/2;
   
   if lmcosi(a1,1) ~= lmcosi(a2,1)
       error('Wrong indeces');
   end
   
   lambda=2*pi/lmcosi(a1,1);
   lambda_linear=lambda*Rref;
   k(i)=1/lambda_linear;   
   
   sdl(i)=2*pi*Rref*Rref*Rref*sum((lmcosi(a1:a2,3).^2)+(lmcosi(a1:a2,4).^2));
        
end 