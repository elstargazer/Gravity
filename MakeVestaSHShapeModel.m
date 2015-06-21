% ccc

tic 

ShapeFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/HASTA_LA_VESTA_SHAPE/SHAPE.TXT';
FiStep=.1;
LambdaStep=.1;

[ri,fii,lambdai]=Load3ColumnsXYZShape(ShapeFileName,FiStep,LambdaStep);

ri(:,end)=ri(:,1);
ri(1,:)=mean(ri(1,:));
ri(end,:)=mean(ri(end,:));

tic
MaxDegree=500;
[lmcosi_shape_s500,dw]=xyz2plm(ri,MaxDegree,'im',[],[],[]);
toc

% Resolution=.05;
% 
% vesta_shape_fig=figure('Position',[1 1 1500 1500],'Color',[0 0 0]);
% plotplm(lmcosi_shape,[],[],2,Resolution,[],[],[]);
%  
toc