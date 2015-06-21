ccc
% d=[2 100:100:1000];

% lmcosi_z=lmcosi_d;
% lmcosi_z(:,3:4)=0;
% lmcosi_z(:,3:4)=0;
% 
% lmcosi_d_tr=lmcosi_d;
% lmcosi_g_tr=lmcosi_g;
% 
% lmcosi_d_tr(:,3:4)=0;
% lmcosi_g_tr(:,3:4)=0;
% 
% 
% progressbar(0);

% Resolution=0.1

% fig=figure('Color','w','Position',[1 1 1400 1400]);
% % hold on;
% 
% 
%  ax=axesm('mollweid','frame','off','FontSize',10,'Grid','on','MLabelParallel',...
%      'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
%      , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',1,'FontColor',[0 0 0]...
%      ,'FontSize',10,'GAltitude',Inf,'Geoid',[1 0]);
% 
% dr_cell=cell(1);
%  
%  
% 
% for i=1:numel(d-1)
%     
%     d1=d(i)+1;
%     d2=d(i+1);
%     
%     s1=d1*(d1+1)/2+1;
%     s2=(d2+1)*(d2+2)/2;
%     
%     lmcosi_d_tr=lmcosi_z;
%     lmcosi_g_tr=lmcosi_z;
%     
%     
%     lmcosi_d_tr(s1:s2,:)=lmcosi_d(s1:s2,:);
%     lmcosi_g_tr(s1:s2,:)=lmcosi_g(s1:s2,:);
%     
%     lmcosi_d_tr(s2+1:end,:)=[];
%     lmcosi_g_tr(s2+1:end,:)=[];
%     
%     r_d_tr=plm2xyz(lmcosi_d_tr);
%     r_g_tr=plm2xyz(lmcosi_g_tr);
%     
%     dr=1000*r_g_tr-r_d_tr;
%     dr_cell{i}=dr;
%     
% %     close all;    
% %     MapRadialGrid(flipud(dr));
%     
%     s=size(dr);
% 
%     fi=linspace(-90,90,s(1))/180*pi;
%     lambda=linspace(0,360,s(2))/180*pi;
%     [lambdai,fii]=meshgrid(lambda,fi);
%     
%     h=surfm(fii,lambdai,flipud(dr));
%     shading interp          
%     
%     title([num2str(d1) '--' num2str(d2)],'FontSize',10);    
%     caxis([-100 100]);
%     set(gca,'FontSize',10);
%     print(gcf,'-depsc',['GaskellDLRdiff' num2str(d1) '-' num2str(d2) '.eps']);
%      
%     unplot
%   
%     progressbar(i/numel(d-1));
%     
% end
% 
% 
% 
% 
% 
% for i=1:numel(dr_cell)
%     
%     d1=d(i)+1;
%     d2=d(i+1);
%     
%     
% %     dr=1000*r_g_tr-r_d_tr;
%     dr=dr_cell{i};
%     
%     
%     s=size(dr);
% 
%     fi=linspace(-90,90,s(1))/180*pi;
%     lambda=linspace(0,360,s(2))/180*pi;
%     [lambdai,fii]=meshgrid(lambda,fi);
%     
%     close all;
%     
%     fig=figure('Color','w','Position',[1 1 1400 1000]);
% 
%     ax=axesm('mollweid','Origin',[0 pi],'frame','off','FontSize',10,'Grid','on','MLabelParallel',...
%      'equator','AngleUnits','radians','LabelUnits','degrees','ParallelLabel'...
%      , 'on','MeridianLabel', 'on','GLineStyle','w--','GlineWidth',1,'FontColor',[0 0 0]...
%      ,'FontSize',10,'GAltitude',Inf,'Geoid',[1 0]);
%  
% %     fii=circshift(fii,[0 fix(s(2)/4)]);
% %     lambdai=circshift(lambdai,[0 fix(s(2)/4)]);
% %     dr=circshift(dr,[0 fix(s(2)/4)]);
%     
% 
%      h=surfm(fii(:,i:end-i),lambdai(:,i:end-i),flipud(dr(:,i:end-i)));
% %     h=surfm(fii,lambdai,flipud(dr));
%     shading interp    
%     
%     
%         
%     cbar=colorbar;
%     
%     ylabel(cbar,'[m]','FontSize',12);
%     
%     title([num2str(d1) '--' num2str(d2)],'FontSize',10);    
%     caxis(0.1*[min(dr(:)) max(dr(:))]);
%     set(gca,'FontSize',10);
%     
%     
%     
%     title([num2str(d1) '--' num2str(d2)],'FontSize',10);    
%     
%     xlabel(['min = ' num2str(min(min(dr(:,i:end-i)))) '; max = ' num2str(max(max(dr(:,i:end-i)))) ' [m]']);
%        
%     print(gcf,'-depsc',['GaskellDLRdiff' num2str(d1) '-' num2str(d2) '.eps']);
%     
%     
% % end


% r_d=plm2xyz(lmcosi_d);
% r_g=plm2xyz(lmcosi_g);


% dr=1000*r_g-r_d;

% MapRadialGrid(flipud(dr))


load Shape1800


[sdl_d,l_d,bta_d,lfit_d,logy_d,logpm_d]=plm2spec(lmcosi_d,2);
[sdl_g,l_g,bta_g,lfit_g,logy_g,logpm_g]=plm2spec(lmcosi_g,2);

% lmcosi_d=lmcosi_g;
lmcosi_diff(:,3:4)=1000*lmcosi_g(:,3:4)-lmcosi_d(:,3:4);

[sdl_diff,l_diff,bta_diff,lfit_diff,logy_diff,logpm_diff]=plm2spec(lmcosi_diff,2);


lambda_wave=2*pi./l_d*265000;

figure; hold on;
set(gca,'FontSize',20);
% set(gca,'XScale','log','Yscale','log')

set(gca,'Yscale','log')

plot(l_d,sdl_d,'b-');
plot(l_d,1000*sdl_g,'r-');

plot(l_d,sdl_diff,'k-');

legend({'DLR','Gaskell','Difference'});

% set(gca,'XDir','reverse')

xlabel('Wavelength [m]','FontSize',20);
ylabel('Spectral density','FontSize',20);


progressbar(1);

%% Computing Correlation

Corr=SphericalHarmonicCorrelation(lmcosi_d,lmcosi_g);

figure; hold on;
set(gca,'FontSize',20);
plot(Corr,'r-','LineWidth',3);
% set(gca,'Yscale','log')


circle_rad=40;
L=15;
MinConcentration=0.1;
Resolution=.5;
MaxDeg=400;

[G2,V2,N2,J2]=glmalphapto(circle_rad,L,0,0);

s=size(G2);

for i=1:s(2)       
    lmcosi_window_basic{i}=glm2lmcosi(G2,i);    
end

Tapers=find(V2>MinConcentration);
NumberOfTapers=numel(Tapers);
lmcosi_window_basic=lmcosi_window_basic(Tapers);

NTapers=sum(V2>MinConcentration);


r_d=plm2xyz(lmcosi_d,Resolution);
r_g=plm2xyz(lmcosi_g,Resolution);

Corr_w=zeros(NumberOfTapers,MaxDeg);

progressbar(0);

for i=1:NumberOfTapers
    
%         [lmcosi_window{i},~,~]=plm2rot(lmcosi_window_basic{i},alp,bta,gam,'dlmb');
    
        [win(:,:,i),~,~,Plm]=plm2xyz(lmcosi_window_basic{i},Resolution);
        r_d_win=r_d.*win(:,:,i);
        r_g_win=r_g.*win(:,:,i);
        
        lmcosi_d_w=xyz2plm(r_d_win,MaxDeg);
        lmcosi_g_w=xyz2plm(r_g_win,MaxDeg);
        Corr_w(i,:)=SphericalHarmonicCorrelation(lmcosi_d_w,lmcosi_g_w);
        
%         size(Corr_w)
        
        progressbar(i/NumberOfTapers);

        
end

progressbar(1);

figure; hold on;
set(gca,'FontSize',20);
ylabel('Correlation []','FontSize',20);

plot(L:MaxDeg-L,Corr_w(end-1,L:end-L)','k-','LineWidth',1);


% MapRadialGrid(win(:,:,4));

































