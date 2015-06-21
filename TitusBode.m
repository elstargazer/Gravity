figure; hold on;
box on
set(gca,'FontSize',20);


a=[0.39 0.72 1.00 1.52 2.77 5.20 9.54 19.23 30.06 39.44];
at = 0.4 + 0.3 * 2 .^ ([-Inf 0:8])

plot([1 2 3 4 5 6 7 8 9 10],a,'x-r','MarkerSize',20)
plot(1:numel(at),at,'o-b','MarkerSize',20)

xlim([1 10])

ylabel('Semimajor axis [AU]')

set(gca,'YScale','log');

ylim([0.2 	78])

set(gca,'XTick',1:10)

xlabel('Object number','FontSize',20);