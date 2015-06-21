function [rinew]=ZonalDifference(ri,riref)

dr=riref-ri;

drnewlonav=mean(dr,2);
% drnewstd=std(drnew,0,2);


rinew=ri;
for i=1:size(rinew,2)
rinew(:,i)=rinew(:,i)+drnewlonav;
end
