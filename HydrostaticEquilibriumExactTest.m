ccc

r=261000;
rho=3457;

Rr=17536000;
Mr=1;


T0=2;
T1=10;
Tstep=0.05;

T=T1:-Tstep:T0;

f0=0;


tic
for i=1:numel(T)       
    [fh(i),d2(i)]=HydrostaticStateExact(r,T(i),rho,f0);
    f0=fh(i);    
end
toc

tic
for i=1:numel(T)       
    [fhr(i),d2r(i)]=HydrostaticStateExactRing(r,T(i),rho,Rr,Mr,f0);
    f0=fhr(i);    
end
toc



figure; hold on;
plot(T,fh,'b-');
plot(T,fhr,'r-');












