function f=RadialGradientFactor(Degree,DerOrder,r)

num=ones(numel(Degree),DerOrder);

for i=1:numel(Degree+1:Degree+DerOrder)
    
    num(:,i)=Degree+i;
    
end

if (DerOrder~=1)
num=(prod(num,2));
end

f=(-1^DerOrder).*num./(r^DerOrder);