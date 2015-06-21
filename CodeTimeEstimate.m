d=[10 20 15 25 30 60 100 150];
t=[4.7 12.88 8.2 19.7 27.9 104.14 356.9 1558.3]


figure; hold on
plot(d,t,'o')

p=polyfit(d,t,3);
d2=2:500;
t2=polyval(p,d2);
% 
% set(gca,'XScale','log')
% set(gca,'YScale','log')

plot(d2,t2,'-b')