ccc


drho=100:100:200;
h=30:5:35;

[drhoi,hi]=meshgrid(drho,h);


bouguer_anomaly=zeros(size(drho));

for i=1:numel(drho)
    for j=1:numel(h)       
        
        
        drhoi(i,j)
        hi(i,j)
        
        [angle_gc_iv,bouguer_anomaly(i,j,:)]=BouguerFromHydrocode(drhoi(i,j),hi(i,j));
        
    end
end