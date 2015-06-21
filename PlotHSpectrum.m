function Hpower=PlotHSpectrum
FileName='file_graph_degvar_total_power';

in=fopen(FileName,'r');

fgets(in);
str=fgets(in);

Hmax=str2num(str(13:end));


for i=1:Hmax
    str=fgets(in);
    Hpower(i)=str2num(str(15:end));
    
    
end


fclose(in);

figure1=PlotHPowers(Hpower);

FileName='hpower.eps';

set(gcf,'InvertHardcopy','off');
print(figure1,'-depsc',FileName)

% set(gcf,'InvertHardcopy','off');
% print(gcf,'-depsc','topopower2.eps')
% 
% set(gcf,'InvertHardcopy','off');
% print(gcf,'-depsc','gravitypower.eps')