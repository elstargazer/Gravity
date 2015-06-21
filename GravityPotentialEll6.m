function [U,lambda]=GravityPotentialEll6(mu,Re,RefEll,lmcosi,x,y,z)

a=RefEll(1);
b=RefEll(2);
c=RefEll(3);

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

s=size(x);

x=x(:);
y=y(:);
z=z(:);

Npts=numel(x);
maxn=lmcosi(end,1);
U=zeros(Npts,1);
Ncoeff=size(lmcosi,1);

%% Converting cartesian to ellipsoidal coordinates
lambda=zeros(Npts,3);

gamma = computeRomainNormalizationConstants(maxn, a, b, c);

lambda = approxCartToEllnosign(a, b, c, [x y z]);    
range_lambda=max(lambda(:,1))-min(lambda(:,1));

[PsiKEven,PsiLEven,PsiMEven,PsiNEven,...
    PsiKOdd,PsiLOdd,PsiMOdd,PsiNOdd]=calcPsiSingle(x,y,z,h,k);

MaxTol=1e-8;

%% Computing potential
progressbar('Computing potential');
                   
for n=0:maxn
            
    ind=n*n+1:(n+1)^2;        
%     alpha=lmcosi(ind,3)';
    
    alpha=repmat(lmcosi(ind,3)'./sqrt(gamma(ind)),Npts,1);  
    
    w=NEllHarm(n);
    
    if (mod(n,2)==0)
        Psi=[repmat(PsiKEven,1,w(1)) repmat(PsiLEven,1,w(2))...
        repmat(PsiMEven,1,w(3)) repmat(PsiNEven,1,w(4))];
    else 
        Psi=[repmat(PsiKOdd,1,w(1)) repmat(PsiLOdd,1,w(2))...
        repmat(PsiMOdd,1,w(3)) repmat(PsiNOdd,1,w(4))];    
    end

    
    [VK, VL, VM, VN] = calcEigCoef(n, h, k);
    [P1,~] = calcPKLMN(lambda(:,1), h, VK, VL, VM, VN);
    [P2,~] = calcPKLMN(lambda(:,2), h, VK, VL, VM, VN);
    [P3,~] = calcPKLMN(lambda(:,3), h, VK, VL, VM, VN);
    
    [E1, Ederiv1] = calcE(lambda(:,2), h, k, VK, VL, VM, VN);
%     [lame,lameDeriv,lameMatrix,lameMatrixDeriv] = calcLame(lambda(:,2), n, 1, a, b, c, 1, 1)
    
    I=zeros(Npts,2*n+1);
    I2=zeros(Npts,2*n+1);
%     Ideriv=zeros(Npts,2*n+1);

  
    if  (range_lambda < MaxTol)
        [I1,~] = computeExteriorIntegralDeg(lambda(1), n, a, b, c);
        for p=1:2*n+1;
            I(:,p)=I1(p);
%             Ideriv(:,p)=Ideriv1(p);
        end
    else
        parfor i=1:Npts
            if isnan(lambda(i))
                I(i,:)=nan(1,2*n+1);
%               Ideriv(i,:)=nan(1,2*n+1);
            else
%                 [I(i,:),~] = calcI(lambda(i,1), n, h, k, VK, VL, VM, VN);
                [I(i,:),~] = computeExteriorIntegralDeg(lambda(i,1), n, a, b, c);
            end
        end
    end
    
    V=Psi.*P1.*P2.*P3.*I;
    
   [Fref,~] = calcLameSecondKind(Re, n, a, b, c, 1, 1);    
   Fref=repmat(Fref,Npts,1);       
  
   U=U+sum(alpha.*V./Fref,2);              
   progressbar((n+1)/(maxn+1));
end 


progressbar(1);
NNAN=sum(isnan(lambda(:,1)));
disp([' Could not find ell. coords in ' num2str(NNAN) ...
    ' or ' num2str(NNAN/Npts*100) '% cases']);

U=reshape(sqrt(8)*mu*U,s);

