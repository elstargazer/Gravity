function [CovMatrix,CSness,Degree,Order]=CovarianceReadErwan(FileName)


% COVMATgeodyn=[COVMATdir,'/',IDfield];

%% open files

fidin=fopen(FileName,'rb');

% read GEODYN-format COVMAT

% HEADER
ntmp=fread(fidin,1,'integer*4','ieee-le');
vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
ntmp=fread(fidin,1,'integer*4','ieee-le');

nparm=round(vtmp(3));
nobs=round(vtmp(8));

% LABELs
ntmp=fread(fidin,1,'integer*4','ieee-le');
vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
ntmp=fread(fidin,1,'integer*4','ieee-le');

lbls=vtmp(1+(1:nparm));

% VALUEs
ntmp=fread(fidin,1,'integer*4','ieee-le');
vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
ntmp=fread(fidin,1,'integer*4','ieee-le');

coefvalues_=vtmp(1+(1:nparm));

% SIGMAs
ntmp=fread(fidin,1,'integer*4','ieee-le');
vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
ntmp=fread(fidin,1,'integer*4','ieee-le');

coefsigmas_=vtmp(1+(1:nparm));

% covariance

covmat_=zeros(nparm,nparm);

for irow=1:nparm
    
    ntmp=fread(fidin,1,'integer*4','ieee-le');
    vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
    ntmp=fread(fidin,1,'integer*4','ieee-le');
    covmat_(irow,irow:nparm)=vtmp(2:end);
    
end

%
fclose(fidin);

%%

lbls=num2str(lbls);

Condition=(lbls(:,1)=='6');

coefvalues=coefvalues_(Condition);
coefsigmas=coefsigmas_(Condition);

CovMatrix=covmat_(:,Condition);
CovMatrix=CovMatrix(Condition,:);

Degree=str2num(lbls(:,7:10)); Degree=Degree(Condition);
Order=str2num(lbls(:,11:14)); Order=Order(Condition);

MaxDegree=max(Degree);
CSness=(Condition & (lbls(:,6)=='1')); 
CSness=CSness(Condition);


CovMatrix(isnan(CovMatrix))=0;

Variance=diag(CovMatrix);
StandardDeviation=sqrt(Variance);




covmat_=covmat_+covmat_'-diag(diag(covmat_));

A=ones(size(covmat_))*diag(diag(covmat_));
corrmat_=covmat_./sqrt(A)./sqrt(A');

% % Make NaN diag
iu = find(triu(ones(nparm,nparm), 0));
corrmat_(iu)=NaN;
    

% corrmat_=triu(corrmat_,1);


[corrsort,ix]=sort(corrmat_(:));


col=fix((ix-1)/nparm)+1;
row=ix-(col-1)*nparm;

% [col row];

corrsort2=corrsort;
for i=1:numel(row)
    corrsort2(i)=corrmat_(row(i),col(i));    
end


in=fopen('CorrelationPars.txt','w');

progressbar(0);

for i=1:numel(row)
    if ~isnan(corrmat_(ix(i)))
        fprintf(in,'%s %s %18.16f\n',lbls(row(i),:),lbls(col(i),:),corrsort(i));
    end    
    progressbar(i/numel(row) );
end

progressbar(1)

fclose(in);

%% Check
clc
i=100000

lbls(row(i),:)
lbls(col(i),:)
corrsort2(i)


a1=lbls(col(i),:);
a2=lbls(row(i),:);

a1='50120070062300';
a2='50120070013200';


i1=1;
i2=1;


for i=1:nparm
    


    if lbls(i,:)==a1
        b1(i1)=i; i1=i1+1;
%         i
    end
    
    if lbls(i,:)==a2
        b2(i2)=i; i2=i2+1;
%         i
    end
    
end

% [b1 b2]


corrmat_(b1(:),b2(:)) 




