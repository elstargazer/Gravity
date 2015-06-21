

ro_mantle=2800:100:3100;
ro_crust=2800:100:3100;


[ro_mantlei,ro_crusti]=meshgrid(ro_mantle,ro_crust);


for i=1:numel(ro_mantle)
    for j=1:numel(ro_crust)
        if (ro_crusti(i,j)<(ro_mantle(i,j)))
            
            [angle_gc,bouguer_anomaly_obs,h_crust]=BouguerAnomalyDH(ro_crust(i,j),ro_mantle(i,j));            
            [angle_gc_iv,bouguer_anomaly_hyd]=BouguerFromHydrocode(delta_rho,h_crust);
            
            R(i,j)=sum((bouguer_anomaly_obs-bouguer_anomaly_hyd).^2);
            
            
            
            
        end
    end
end