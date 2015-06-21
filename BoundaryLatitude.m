function fi=BoundaryLatitude(lambda,lattrk,lontrk)

[lontrk,m,n]=unique(lontrk);


lattrk=lattrk(m);


fi=interp1(lontrk,lattrk,lambda,'nearest');

