function [binMean,binStandardDev]=2Drebin( x,y,z,sx,sy )


mask=~(isnan(x) | isnan(y) | isnan(z));
x=x(mask);
y=y(mask);
z=z(mask)

XbotEdge=min(x(:));
XtopEdge=max(x(:));

YbotEdge=min(x(:));
YtopEdge=max(x(:));

XbinEdges=XbotEdge:sx:XtopEdge;
YbinEdges=YbotEdge:sy:YtopEdge;

[~,XwhichBin]=histc(x,XbinEdges);
[~,YwhichBin]=histc(y,YbinEdges);

% progressbar(0)

for i=0:numel(XbinEdges)
    for j=0:numel(YbinEdges)
    
    XflagBinMembers=(XwhichBin==i);
    YflagBinMembers=(YwhichBin==i);
    
    binMembers=z(XflagBinMembers & YflagBinMembers);
    
    binMean(i+1,j+1)=mean(binMembers);  
    binStandardDev(i+1,j+1)=std(binMembers);  
    
%     progressbar(i/numel(binEdges));
    
    end   
end
% 
% progressbar(1)
% 
% binMean=binMean(2:end);
% binStandardDev=binStandardDev(2:end);
% binCenters=binEdges+s/2;




