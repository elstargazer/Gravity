function [XbinCenters,YbinCenters,binMean,binStandardDev]=rebin2D( x,y,z,z_std,sx,sy )


mask=~(isnan(x) | isnan(y) | isnan(z));
x=x(mask);
y=y(mask);
z=z(mask);

XbotEdge=min(x(:));
XtopEdge=max(x(:));

YbotEdge=min(y(:));
YtopEdge=max(y(:));

XbinEdges=XbotEdge:sx:XtopEdge;
YbinEdges=YbotEdge:sy:YtopEdge;

[~,XwhichBin]=histc(x,XbinEdges);
[~,YwhichBin]=histc(y,YbinEdges);

% progressbar(0)

for i=0:numel(XbinEdges)
    for j=0:numel(YbinEdges)
    
    XflagBinMembers=(XwhichBin==i);
    YflagBinMembers=(YwhichBin==j);
    
    binMembers=z(XflagBinMembers & YflagBinMembers);
    binMembers_std=z_std(XflagBinMembers & YflagBinMembers);
    
    binMembers_weight=(1./binMembers_std).^2;
    

    %binMean(i+1,j+1)=mean(binMembers);
    
    binMean(i+1,j+1)=sum(binMembers.*binMembers_weight)./sum(binMembers_weight);
       
    binStandardDev(i+1,j+1)=std(binMembers);  
    
%     progressbar(i/numel(binEdges));
    
    end   
end
% 
% progressbar(1)
% 
binMean=binMean(2:end,2:end);
binStandardDev=binStandardDev(2:end,2:end);

XbinCenters=XbinEdges+sx/2;
YbinCenters=YbinEdges+sy/2;




