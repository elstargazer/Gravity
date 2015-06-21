function [dr_db,fi_db]=DebiasHist(dr,fi)


n=size(dr,1);


cfi=cos(fi(:,1)/180*pi);

npts=fix(cfi*n);

dr_db=nan(1,numel(dr));
fi_db=dr_db;

ind_start=1;

progressbar(0);


for i=2:n-1
    
    ind=fix(rand(1,npts(i))*n)+1;
    
    ind_finish=ind_start+numel(ind)-1;
    
%     [ind_start ind_finish]
%     i
%     
    dr_db(ind_start:ind_finish)=dr(i,ind);
    fi_db(ind_start:ind_finish)=fi(i,ind);
    
%     [ind_start ind_finish]
    
    ind_start=ind_finish+1;
    
    progressbar(i/(n-2));
    
    
    
end

progressbar(1);

dr_db(isnan(dr_db))=[];
fi_db(isnan(fi_db))=[];