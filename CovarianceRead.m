function [CovMatrix,CSness,mu,lmcosi,Rref,Degree,Order]=CovarianceRead(FileName)
%  ccc

%% Covariance Matrix File
% FileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.DAT';
% FileName='/Users/antonermakov/Moon/gglp_glgm3_shb.dat';


in=fopen(FileName,'r','l','US-ASCII');


Rref=fread(in,1,'float64');
mu=fread(in,1,'float64');
mu_std=fread(in,1,'float64');
MaxDegree=fread(in,1,'int32');
MaxOrder=fread(in,1,'int32');
NormalizationState=fread(in,1,'int32');
NumberOfNames=fread(in,1,'int32');
RefLon=fread(in,1,'float64');
RefLat=fread(in,1,'float64');
 


%% Find the black area
NumberOfSpaces=57;
Blank=fread(in,8*NumberOfSpaces,'*char');


CoefName=fread(in,[8 NumberOfNames],'*char')';


%% Find the black area
NumberOfSpaces=10;
Blank=fread(in,8*NumberOfSpaces,'*char');

Coef=fread(in,NumberOfNames,'float64');


CovSize=0.5*NumberOfNames*(NumberOfNames+1);

%% Find the black area
NumberOfSpaces=10;
Blank=fread(in,8*NumberOfSpaces,'*char');

Cov=fread(in,CovSize,'float64');


CovMatrix=zeros([NumberOfNames NumberOfNames]);
% CorrMatrix=zeros([NumberOfNames NumberOfNames]);

progressbar('1/4 step: covariance matrix');

%% Covariance Matrix

for i=1:NumberOfNames
    
    CovMatrix(1:i,i)=Cov(0.5*i*(i-1)+1:0.5*i*(i+1));   
    
%     start=0.5*i*(i-1)+1
%     finish=0.5*i*(i+1)  

progressbar(i/NumberOfNames);

end

Variance=diag(CovMatrix);
StandardDeviation=sqrt(Variance);

clear Cov

progressbar('2/4 step: correlation matrix');

%% Correlaton Matrix


% for i=1:NumberOfNames
%    for j=1:i
%        CorrMatrix(j,i)=CovMatrix(j,i)/(StandardDeviation(i)*StandardDeviation(j));
%    end
%    
%    progressbar(i/NumberOfNames);
% end



progressbar('3/4 step: degrees and orders')

%% Degrees and orders 

Degree=str2num(CoefName(2:end,2:4));
Order=str2num(CoefName(2:end,5:7));
Type=('C'==(CoefName(2:end,1)));

Coef(1)=[];

[~,~,~,lmcosi,~,~,~,~,~,~,~]=addmon(MaxDegree);

progressbar('4/4 step: model coefficients')

%% Model coefficients

for i=1:numel(Coef)
    
    
    row=0.5*Degree(i)*(Degree(i)+1)+1+Order(i);    
    if (Type(i)==1)        
        lmcosi(row,3)=Coef(i);
    else
        lmcosi(row,4)=Coef(i);        
    end
    
    progressbar(i/numel(Coef));
end
lmcosi(1,3)=mu;
fclose(in);


%% Proceesing Cov Matrix

CSness=CoefName(2:end,1)=='C';
CovMatrix(1,:)=[];
CovMatrix(:,1)=[];


