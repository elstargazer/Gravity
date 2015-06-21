ccc;
% 
filename='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/Corrected/SHAPE_NEW.TXT';
FiStep=1.8;
LambdaStep=1.8;

[ri_gas,fii_gas,lambdai_gas]=Load3ColumnsXYZShape(filename,FiStep,LambdaStep);

ri_gas=ri_gas(:);
fii_gas=fii_gas(:);
lambdai_gas=lambdai_gas(:);

MaxDegree=100;
% tic;
% [lmcosi,dw]=xyz2plm(r2,MaxDegree,'irr',fi2*180/pi,lambda2*180/pi,1e-2);
% toc
tic;
[lmcosi,dw]=xyz2plm(ri_gas,MaxDegree,'im',fii_gas*180/pi,lambdai_gas*180/pi,[]);
toc;
