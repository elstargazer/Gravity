function U=GravityPotentialTaylor(lmcosi,Rref,mu,r_b,MaxDerOrder,x,y,z)

[Lambda,Fi,r]=cart2sph(x,y,z);

[x_b,y_b,z_b]=sph2cart(Lambda,Fi,r_b);

U_b=GravityPotential(mu,Rref,lmcosi,x_b,y_b,z_b);

U_taylor=U_b.*0;

for DerOrder=1:MaxDerOrder
    
lmcosi_grad=ComputeGradientCoefficients(lmcosi,DerOrder,r_b);

U_der=GravityPotential(mu,Rref,lmcosi_grad,x_b,y_b,z_b);

mean(mean(U_der));

U_taylor=U_taylor+U_der.*((r-r_b).^DerOrder)./factorial(DerOrder);

% max(max(abs(U_der.*((r-r_b).^DerOrder)./factorial(DerOrder))))


% surf(x,y,z,U_b+U_taylor);
% shading interp
% lighting phong
% axis equal
% colorbar
% drawnow



end

U=U_b+U_taylor;