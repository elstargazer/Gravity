function d2=DeltaSqGravity(rho,xi,yi,zi,FV,xr,yr,zr,gx,gy,gz)

[ax,ay,az]=GravityAccelerationTriDen(xi',yi',zi',FV(:,1)',FV(:,2)',FV(:,3)',xr,yr,zr,rho);

d2=sum((ax-gx).^2+(ay-gy).^2+(az-gz).^2);