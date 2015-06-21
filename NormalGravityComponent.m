function [gn]=NormalGravityComponent(gx,gy,gz,x,y,z)

g=sqrt(gx.*gx+gy.*gy+gz.*gz);

gx_u=gx./g;
gy_u=gy./g;
gz_u=gz./g;


[Nx,Ny,Nz] = surfnorm(x,y,z);


N=sqrt(Nx.*Nx+Ny.*Ny+Nz.*Nz);

Nx_u=-Nx./N;
Ny_u=-Ny./N;
Nz_u=-Nz./N;

l_n=gx_u.*Nx_u+gy_u.*Ny_u+gz_u.*Nz_u;

gn=l_n.*g;