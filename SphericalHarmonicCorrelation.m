function Corr=SphericalHarmonicCorrelation(lmcosi1,lmcosi2)

MaxDegree=min(lmcosi1(end,1),lmcosi2(end,1));

Corr=zeros(1,MaxDegree);

for i=1:MaxDegree+1
    
    i_start=(i-1)*(i)/2+1;
    i_end=(i)*(i+1)/2;   
        
    Corr(i)=sum(lmcosi1(i_start:i_end,3).*lmcosi2(i_start:i_end,3)+lmcosi1(i_start:i_end,4).*lmcosi2(i_start:i_end,4))./...
        sqrt(sum(lmcosi1(i_start:i_end,3).^2+lmcosi1(i_start:i_end,4).^2)*sum(lmcosi2(i_start:i_end,3).^2+lmcosi2(i_start:i_end,4).^2));
end

% figure;
% plot(0:MaxDegree,Corr,marker,'Color',color,'LineWidth',2,'MarkerSize',2);
% title('Gravity/Topography Correlation','FontSize',25,'FontName','Times');
% xlabel('Degree','FontSize',25,'FontName','Times');

% grid on;hold on;

