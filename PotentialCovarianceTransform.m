function [CovMatrixT,der]=PotentialCovarianceTransform(FileName,x,y,z,TruncationDegree)

%% load gravity model and covariances

[CovMatrix,CSness,mu,lmcosi,Rref,Degree,Order]=CovarianceRead(FileName);

% [lmcosi,mu,Rref]=ReadGEODYNGravity('/Users/antonermakov/Dawn/Gravity/Frank20/grvfld.dawn.img.test37');

%% Truncating

TruncationCondition=Degree>TruncationDegree;

Degree(TruncationCondition)=[];
Order(TruncationCondition)=[];
CovMatrix(TruncationCondition,:)=[];
CovMatrix(:,TruncationCondition)=[];
CSness(TruncationCondition)=[];
lmcosi=TruncateGravityModel(lmcosi,TruncationDegree,1);

Ncoeff=numel(CSness);
Npts=numel(x);

%% Deleting first row and column of cov. matrix => GM covariances

CovMatrix=CovMatrix+CovMatrix'-diag(diag(CovMatrix));

% r=sqrt(x.*x+y.*y+z.*z);

mu=mu*1e9;
lmcosi(1,3)=1;
Rref=Rref*1000;
d2=zeros(Ncoeff,Npts);

parfor i=1:Npts
    
    if (isnan(x(i)))        
        d2(:,i)=NaN;        
    else        
        [d,CSness2,Degree2,Order2]=GravityPotentialDerivatives(mu,Rref,lmcosi,x(i),y(i),z(i));

        [~, P] = ismember([CSness Degree Order],[CSness2 Degree2 Order2],'rows');
        der(:,i)=d(P);
    end
end

CovMatrixT=CovarianceTransform(CovMatrix,der);
