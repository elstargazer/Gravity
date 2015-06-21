ccc

a=3;
b=2;
c=1;


k=sqrt(a*a-c*c);
h=sqrt(a*a-b*b);

lambda=-k:0.01:k;


signm=sign(abs(lambda)-h);
signn=sign(abs(lambda)-k);

n=3;


parfor i=1:numel(lambda)
[En(i,:),~] = calcN(lambda(i), n, a, b, c,signm(i), signn(i));
end

Eplotn=En;

n=3;p=7;

Ea=LameFirstKindAn(lambda',h,k);
Eplota=(Ea(:,EllIndex(n,p)));


figure; hold on;

plot(lambda/k,Eplotn,'r');
plot(lambda/k,Eplota,'b');


% figure; hold on;
% plot(diff(Eplotn),'r');


% figure; hold on;
% plot(Eplotn,Eplota,'b');

xlim([-1 1]);

% n=10;
% p=1;
% 
% lambda=a;
% mu=h^2;
% nu=0;
% 
% [lame,lameDeriv,lameMatrix,lameMatrixDeriv] = calcLame(a, n, p, a, b, c, ...
% 						  mu, nu);
%                       
%                       
% size(lame)
% 
% size(lameMatrix)