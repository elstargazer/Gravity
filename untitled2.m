
% 
% tol1=.1;
% eps1=2;
% 
% tol2=0.01;
% eps2=2;
% 
% 
% radii=[490000; 490000; 460000];
% cg_req=max(r_grid(:));
% 
% [xg,yg,zg]=TriEllRadVec(fi_grid',lambda_grid',radii(1),radii(1),radii(3),'xyz');
% U0=GravityPotential(mu,Rref_grav,lmcosi_grav,xg,yg,zg);
% while eps2>tol2 
%     eps1=2;
%     while eps1>tol1
%         
%         [xg,yg,zg]=sph2cart(lambda_grid',fi_grid',rg);
%         
%         U=GravityPotential(mu,Rref_grav,lmcosi_grav,xg,yg,zg);
%         rg=sqrt(xg.*xg+yg.*yg+zg.*zg);
%         
%         Urot=0.5*(omega.^2)*(xg.*xg+yg.*yg);
%         
%         U=U+Urot;
%         
%         g=mu/(mean(rs(:)).^2);
%         
%         dr=-(U0-U)/g;
%         
%         rg=rg+dr;
%         
%         eps1=max(abs(dr(:)));
%     end
%     radii_old=radii;
%     
%     [ ~, radii, ~, ~ ] = ellipsoid_fit( [xg(:) yg(:) zg(:)], 4 )
%     
%     eps2=max(abs(radii_old-radii)) 
%     U0=U0*radii(3)/cg_req;
% end
