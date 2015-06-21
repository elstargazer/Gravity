function BalminoVSSimons(lmcosi_grav,lmcosi_grav_b,Rref,mu)

% lmcosi_grav_b=ReadBalminoSH('file_harmo_pot');
ag=2.918299428885863e+05;
bg=2.650067859489697e+05;
[xe,ye,ze]=MakeRotationalEllipsoid(ag,bg,1,1);
% Rref=265000;
[gx,gy,gz]=GravityAcceleration(mu,Rref,lmcosi_grav,xe,ye,ze);
[gx_b,gy_b,gz_b]=GravityAcceleration(mu,Rref,lmcosi_grav_b,xe,ye,ze);
[g_up,g_east,g_north]=GravityComponents(gx,gy,gz,xe,ye,ze,ag,bg);
[g_up_b,g_east_b,g_north_b]=GravityComponents(gx_b,gy_b,gz_b,xe,ye,ze,ag,bg);


a=1;
b=-1;

MapRadialGrid((a*g_up+b*g_up_b)'*1e5); 
