R_R=450/2/265*180/pi;
R_V=395/2/265*180/pi;

fi_R=-75;
lambda_R=301;
fi_V=-52;
lambda_V=170;

% [latc_R,lonc_R] = scircle1(fi_R,lambda_R,R_R);
% [latc_V,lonc_V] = scircle1(fi_V,lambda_V,R_V);
% 
% plotm(latc_R/180*pi,lonc_R/180*pi,'--k','LineWidth',4);
% plotm(latc_V/180*pi,lonc_V/180*pi,'--k','LineWidth',4);
% 

PlotEllipseOnSphere(-19,-130,33,20,90); % PlotVestalia Terra
PlotEllipseOnSphere(fi_R,lambda_R,R_R,R_R,0); % Plot Rheasilvia 
PlotEllipseOnSphere(fi_V,lambda_V,R_V,R_V,0); % Plot Veneneia 