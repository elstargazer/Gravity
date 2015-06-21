function [g_up,g_east,g_north]=GravityComponents(gx,gy,gz,x,y,z,ag,bg)


[B,L,H]=XYZ2BLH(x,y,z,ag,Eccentricity(ag,bg));


North_x=-cos(L).*sin(B);
North_y=-sin(L).*sin(B);
North_z=cos(B);

% quiver3(x,y,z,North_x,North_y,North_z,'r')
% hold on;

% North=sqrt(North_x.*North_x+North_y.*North_y+North_z.*North_z);


East_x=-sin(L);
East_y=cos(L);
East_z=0.*L;

% quiver3(x,y,z,East_x,East_y,East_z,'g')


% East=sqrt(East_x.*East_x+East_y.*East_y+East_z.*East_z);

Up_x=cos(L).*cos(B);
Up_y=sin(L).*cos(B);
Up_z=sin(B);

% Up=sqrt(Up_x.*Up_x+Up_y.*Up_y+Up_z.*Up_z);

% quiver3(x,y,z,Up_x,Up_y,Up_z,'b')

% axis equal

g=sqrt(gx.*gx+gy.*gy+gz.*gz);

gx_u=gx./g;
gy_u=gy./g;
gz_u=gz./g;

l_up=gx_u.*Up_x+gy_u.*Up_y+gz_u.*Up_z;
l_east=gx_u.*East_x+gy_u.*East_y+gz_u.*East_z;
l_north=gx_u.*North_x+gy_u.*North_y+gz_u.*North_z;

g_up=l_up.*g;
g_east=l_east.*g;
g_north=l_north.*g;

% g2=sqrt(g_up.*g_up+g_north.*g_north+g_east.*g_east);
% 
% dg=g2-g;

% min(min(dg))
% max(max(dg))