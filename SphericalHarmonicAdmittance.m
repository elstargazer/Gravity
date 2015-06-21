function Adm=SphericalHarmonicAdmittance(lmcosi1,lmcosi2)

MaxDegree=min(lmcosi1(end,1),lmcosi2(end,1));

Adm=zeros(1,MaxDegree);



for i=3:MaxDegree+1;
    
    i_start=(i-1)*(i)/2+1;
    i_end=(i)*(i+1)/2;   
    
    
    
    Adm(i-1)=sum(lmcosi1(i_start:i_end,3).*lmcosi2(i_start:i_end,3)+lmcosi1(i_start:i_end,4).*lmcosi2(i_start:i_end,4))./...
        sum(lmcosi2(i_start:i_end,3).^2+lmcosi2(i_start:i_end,4).^2);
         
end

% n=1:MaxDegree;


% figure;
% plot(n,Adm,'o-','Color',color,'LineWidth',4,'MarkerSize',2);
% title('Gravity/Topography Admittance','FontSize',25,'FontName','Times');
% xlabel('Degree','FontSize',25,'FontName','Times');
% 
% grid on;
% hold on;
% 
