ccc

x=rand(1,100);


n = 2;
m = length(x);
%% Generate the z variable as a mapping of your x data range into the 
%% interval [-1,1]

z = ((x-min(x))-(max(x)-x))/(max(x)-min(x));

A(:,1) = ones(m,1);
if n > 1
   A(:,2) = z;
end
if n > 2
  for k = 3:n+1
     A(:,k) = 2*z.*A(:,k-1) - A(:,k-2);  %% recurrence relation
  end
end

b = A \ y;
