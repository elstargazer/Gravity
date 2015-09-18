function [param]=mcmc_2l_run(N, s, data, param0, cost_fun,step_param,T,r2,M)

param=zeros(N,2);

param(1,:)=param0;

cost_c=cost_fun(param(1,:),data,T,r2,M);

progressbar(0);

for i=1:N
    
    param_p = step_param(param(i,:),s);
    
    
    r1_p=param_p(1);
    rho1_p=param_p(2);
   
    rho2=-(3*M-4*pi*(r1_p.^3).*rho1_p)./...
    (4*pi*(r1_p.^3)-4*pi*(r2^3));

    cond = ((rho2>0) && (rho1_p>0) && (r1_p>0));
    
    if cond
        cost_p = cost_fun(param_p,data,T,r2,M);
    else
        param_p = param0;
        cost_p = cost_fun(param_p,data,T,r2,M);
        1
    end
    
    metric=exp(- (cost_p) + (cost_c));
     
    if (rand < metric) && ~isnan(cost_p) && ~isnan(metric)
        
        param(i+1,:) = param_p;    
        cost_c=cost_p;
        
    else
        param(i+1,:) = param(i,:);  

    end
    
    progressbar(i/N);
    
end


 progressbar(1);