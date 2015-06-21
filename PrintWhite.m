function PrintWhite(FileName)

set(gcf,'color','w');
set(gcf,'InvertHardcopy','on');
print(gcf,'-depsc',FileName);