function R=MeanEllRadius(a,b,c)

R=zeros(size(a));

for i=1:numel(a)  
    
    fun = @(fi,lambda) cos(fi).*TriEllRadVec(fi,lambda,a(i),b(i),c(i),'rad');

    R(i) = integral2(fun,-pi/2,pi/2,0,2*pi)/(4*pi);  
end

