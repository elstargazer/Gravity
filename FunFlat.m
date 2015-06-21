function df=FunFlat(rhocore,M,rcore,router,T,fouter_obs)

fouter0=0.1;
fcore0=0.1;


rhoouter=-(3*M-4*pi*(rcore^3)*rhocore)/(4*pi*(rcore^3)-4*pi*(router^3));
    
[fh,~]=HydrostaticStateExact2l(router,rcore,T,rhoouter,rhocore,fouter0,fcore0);
            
            
df=fh(1)-fouter_obs;