% ccc
tic

%% Input parameters

%  FileName='/Users/antonermakov/Dawn/Gravity/VESTA20G/JGV20G02.DAT';
CovFileName='/Users/antonermakov/Dawn/Gravity/Frank20/covmat.surv_ham2conantwpim';
GrModelFileName='/Users/antonermakov/Dawn/Gravity/Frank20/grvfld.dawn.img.surv_ham2conantwpim_41';

TruncationDegree=8; % gravity model truncation degree

StepL=10; % contour step in mGal

FiStep=5; % latitude step
LambdaStep=5; % longitude step


ag=290000; % reference ellipsoid semimajor axis
bg=265000; % reference ellipsoid semiminor axis



%% load gravity model and covariances

% [CovMatrix,CSness,mu,lmcosi,Rref,Degree,Order]=CovarianceRead(FileName);

[CovMatrix,CSness,Degree,Order]=CovarianceReadErwan(CovFileName);
[lmcosi,mu,Rref]=ReadGEODYNGravity(GrModelFileName);



% figure; hold on;
% set(gca,'FontSize',20,'YScale','Log');
% xlim([2 TruncationDegree]);
% 
% 
% [sdl,l,~,~,~,~]=plm2spec(lmcosi);
% plot(l,sdl,'-o','LineWidth',3,'MarkerSize',10);
% 
% [sdl,l,~,~,~,~]=plm2spec([lmcosi(:,1:2) lmcosi(:,5:6)]);
% plot(l,sdl,'--ro','LineWidth',3,'MarkerSize',10);
% 
% 
% 
% 
% Sigma=sqrt(diag(CovMatrix));
% 
% nm=[lmcosi(:,1:2); lmcosi(:,1:2)]
% 
% [~, P] = ismember([Degree Order],[],'rows');
% 
% 
% for i=1:size(lmcosi,1)
%     
%     
%     [~, P] = ismember([lmcosi(i,1:2) 1],[Degree Order CSness],'rows'); 
%     if P~=0
%         lmcosi(i,5)=Sigma(P);
%     end
%     
%     [~, P] = ismember([lmcosi(i,1:2) 1],[Degree Order CSness],'rows');    
%     if P~=0
%         lmcosi(i,6)=Sigma(P);   
%     end
%     
% end
% 
% [sdl,l,~,~,~,~]=plm2spec([lmcosi(:,1:2) lmcosi(:,5:6)]);
% plot(l,sdl,'--ko','LineWidth',3,'MarkerSize',10);




%% Truncating

TruncationCondition=Degree>TruncationDegree;

Degree(TruncationCondition)=[];
Order(TruncationCondition)=[];
CovMatrix(TruncationCondition,:)=[];
CovMatrix(:,TruncationCondition)=[];
CSness(TruncationCondition)=[];
lmcosi=TruncateGravityModel(lmcosi,TruncationDegree,1);
Ncoeff=numel(CSness);

%% Deleting first row and column of cov. matrix => GM covariances

CovMatrix=CovMatrix+CovMatrix'-diag(diag(CovMatrix));

A=ones(size(CovMatrix))*diag(diag(CovMatrix));

CorrMatrix=CovMatrix./sqrt(A)./sqrt(A');

% plot correlation matrix

figure; hold on;
set(gca,'FontSize',20);
imagesc(CorrMatrix);
caxis([-1 1]);
cbar=colorbar('FontSize',20);
xlim([1 size(CorrMatrix,1)]);
ylim([1 size(CorrMatrix,2)]);

ylabel('Parameter Number','FontSize',20);
xlabel('Parameter Number','FontSize',20);
ylabel(cbar,'Correlation Coefficient []','FontSize',20);



%% Reference surface

[lambdai,fii]=meshgrid(0:LambdaStep/180*pi:2*pi, -pi/2:FiStep/180*pi:pi/2 );

x=ag*cos(lambdai).*cos(fii);
y=ag*sin(lambdai).*cos(fii);
z=bg*sin(fii);

Npts=numel(x);

%% Reference surface, plane where gravity is computed

% Shape model in SH

% load VestaHASTALAVESTAshape_sh720.mat;
% MaxDegreeTopo=50;
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);

% y=-400000:5000:+400000;
% z=-400000:5000:+400000;
% 
% [y,z]=meshgrid(y,z);
% 
% x=0.*y;
% 
% s_grid=size(x);
% 
% [lambda_surf,fi_surf,r_surf]=cart2sph(x,y,z);
% 
% [r_topo_sh_cond,~,~,~]=plm2xyz(lmcosi_shape,fi_surf(:)*180/pi,lambda_surf(:)*180/pi);
% 
% r_topo_sh_cond=reshape(r_topo_sh_cond,s_grid);
% 
% condition=(r_surf>r_topo_sh_cond);
% 
% 
% x(~condition)=NaN;
% y(~condition)=NaN;
% z(~condition)=NaN;

%% Computing errors in potential and total acceleration

sigma_U=zeros(size(x));
sigma_a_rad=zeros(size(x));
% 
% x=300000;
% y=270000;
% z=320000;

r=sqrt(x.*x+y.*y+z.*z);

% progressbar('Computing errors');


% The following three lines are not required for GEODYN models
mu=mu*1e9;
lmcosi(1,3)=1;
Rref=Rref*1000;

% i=1;

U=GravityPotential(mu,Rref,lmcosi,x,y,z);
[gx_obs,gy_obs,gz_obs]=GravityAcceleration(mu,Rref,lmcosi,x,y,z);

g_obs=sqrt(gx_obs.*gx_obs+gy_obs.*gy_obs+gz_obs.*gz_obs);


tic
parfor i=1:numel(x)
    
    if (isnan(x(i)))        
        sigma_a_rad(i)=NaN;        
    else        
        [d,CSness2,Degree2,Order2]=GravityPotentialDerivatives(mu,Rref,lmcosi,x(i),y(i),z(i));

        [~, P] = ismember([CSness Degree Order],[CSness2 Degree2 Order2],'rows');

        d=d(P);
        d2=d.*(Degree+1)./r(i);
        % sigma_U(i)=(d')*CovMatrix*(d);
        % d2=[g_rad_obs(i)/mu; d2];
        sigma_a_rad(i)=(d')*CovMatrix*(d);

    end
end
t_cov=toc

sigma_a_rad=sqrt(sigma_a_rad);


%% Plotting errors in radial gravity acceleration in mGals

MinL=fix(min(sigma_a_rad(:)));
MaxL=fix(max(sigma_a_rad(:)))+1;


% labels=MinL:StepL:MaxL;

AGUaxes;
pcolorm(fii,lambdai,sigma_a_rad); shading interp;
cbar=colorbar('FontSize',20);
ylabel(cbar,'Acceleration standard deviation [mGal]','FontSize',20);
% [cont_fig,h]=contourm(fii,lambdai,sigma_a_rad,labels,'Color','k','LineWidth',3);
% clabelm(cont_fig,h,labels);

% figure
% pcolor(y,z,log10(sigma_a_rad)); shading interp; colorbar
% save('sigma_a.mat','sigma_a_rad','y','z');

%% Eigensystem


[V,D] = eig(CovMatrix);





toc