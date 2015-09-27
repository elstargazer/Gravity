function lmcosi_filtered = ApplyCosFilter(lmcosi,n1,n2)

lmcosi_filtered = lmcosi;
n = lmcosi(:,1);
f = filter(n,n1,n2);
lmcosi_filtered(:,3:4) = lmcosi_filtered(:,3:4).*[f f];
lmcosi_filtered = TruncateGravityModel(lmcosi_filtered,n2-1,1);


function f = filter(n,n1,n2)

f = zeros(numel(n),1);

for i=1:numel(n)
    
    if n(i)<n1
        f(i) = 1;
    elseif n(i)>n2
        f(i) = 0;
    else
        f(i)=0.5*(1-cos(pi*n(i)/(n2-n1)));
    end
    
end