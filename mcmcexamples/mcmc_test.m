function [param]=mcmc_test(N, s, data, param0, cost_fun,step_param)

param=zeros(N,3);

param(1,:)=param0;

cost_c=cost_fun(param(1,:),data);

for i=1:N
    
    param_p = step_param(param(i,:),s);  
    cost_p = cost_fun(param_p,data);
    
    metric=exp(- log(cost_p) + log(cost_c));
    
    cond=((param_p(1)>param_p(2)) & (param_p(2)>param_p(3)));
    
%     m(i)=metric;
    
    if (rand < metric) && cond
        
        param(i+1,:) = param_p;    
        cost_c=cost_p;  
        
    else
        param(i+1,:) = param(i,:);      
    end
    
    progressbar(i/N)
    
end