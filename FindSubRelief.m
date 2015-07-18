function lmcosi_sub = FindSubRelief(lmcosi_grav,lmcosi_shape,GM,rho1,rho2)

nc = 5;
eps_ri_mantle = 1;

Coef1=M*(2.*n+1).*((Rref./D_mantle).^n)./...
    (4*pi*ro_mantle_diff.*(D_mantle^2));

ll=1./(M*(2.*nc+1).*((Rref./D_mantle).^nc)./...
    (4*pi*ro_mantle_diff.*(D_mantle^2))).^2;

w=(1+ll.*(Coef1.^2)).^(-1);

lmcosi_mantle_shape_correction_1(:,3:4)=[w,w].*(lmcosi_ba(:,3:4)).*[Coef1,Coef1];

lmcosi_mantle_shape_correction_full=lmcosi_mantle_shape_correction_1;

dri_mantle=eps_ri_mantle+1;

iter=1;

while (dri_mantle>eps_ri_mantle)
    
    [ri_mantle_shape_correction,~]=plm2xyz(lmcosi_mantle_shape_correction_full,1);
    
    n=2;
    
    dri_mantle_fa_correction=eps_mantle_fa_correction+1;
    
    while (dri_mantle_fa_correction>eps_mantle_fa_correction)
        
        P=ones(numel(Degree),1);
        
        for j=1:n
            P=P.*(Degree+4-j);
        end
        
        Coef2=P./((D_mantle^n)*factorial(n).*(Degree+3));
        
        ri_mantle_shape_correction_n=ri_mantle_shape_correction.^n;
        
        [lmcosi_mantle_shape_correction_power,dw]=xyz2plm(ri_mantle_shape_correction_n,N_trunc,'im',[],[],[]);
        
        lmcosi_mantle_shape_fa_correction(:,3:4)=[w, w].*D_mantle.*lmcosi_mantle_shape_correction_power(:,3:4).*[Coef2, Coef2];
        
        lmcosi_mantle_shape_correction_full(:,3:4)=lmcosi_mantle_shape_correction_1(:,3:4)-lmcosi_mantle_shape_fa_correction(:,3:4);
        
        % lmcosi_shape_mantle_new(:,3:4)=lmcosi_shape_mantle(:,3:4)+lmcosi_shape_correction(:,3:4);
        %
        % plotplm(lmcosi_shape_mantle_new,[],[],2,1,[],[],[]);
        [ri_mantle_fa_correction,~]=plm2xyz(lmcosi_mantle_shape_fa_correction,1);
        
        dri_mantle_fa_correction=max(max(abs(ri_mantle_fa_correction)));
        
        n=n+1;
        
    end
    
    lmcosi_mantle_shape_new(:,3:4)=lmcosi_mantle_shape(:,3:4)+lmcosi_mantle_shape_correction_full(:,3:4);
    
    % plotplm(lmcosi_mantle_shape_new,[],[],2,1,[],[],[]);
    [ri_mantle_new,~]=plm2xyz(lmcosi_mantle_shape_new,1);
    
    dri_mantle=max(max(abs(ri_mantle_new-ri_mantle)));
    
    ri_mantle=ri_mantle_new;
    
    iter=iter+1;
    
end