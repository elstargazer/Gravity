%%
clear

COVMATdir='/Users/antonermakov/Dawn/Gravity/Frank20';
IDfield='covmat.dawn.img.test37';
COVMATgeodyn=[COVMATdir,'/',IDfield];

%% open files

fidin=fopen(COVMATgeodyn,'rb');

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

covmat_=NaN*ones(nparm,nparm);

for irow=1:nparm
    
    ntmp=fread(fidin,1,'integer*4','ieee-le');
    vtmp=fread(fidin,ntmp/8,'real*8','ieee-le');
    ntmp=fread(fidin,1,'integer*4','ieee-le');
    covmat_(irow,irow:nparm)=vtmp(2:end);
    
end

%
fclose(fidin);