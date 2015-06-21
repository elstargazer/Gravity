function PrintBlack(FileName)
% set(gcf,'PaperType','A3')

% set(gcf, 'Units','centimeters', 'Position',[0 0 40 40])
% set(gcf, 'PaperPositionMode','auto')

set(gcf,'color','k');
set(gcf,'InvertHardcopy','off');
% set(gcf,'Position',[1 1 1200 1200]);
print(gcf,'-depsc',FileName);