T=4.4000

r1=273.8448
r2=252.4160
r3=122.5902

rho1=2858.7207
rho2=3316.5337
rho3=7850.5728

f0=[0.1 0.1 0.1]

[fh,d2]=HydrostaticStateExact3l(r1,r2,r3,T,rho1,rho2,rho3,f0(1),f0(2),f0(3));

fh'