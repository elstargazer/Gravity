function [rinew,delta,angle_min]=FindLongitudeOffset(lambdai,ri,riref,exclude,...
    min_angle,max_angle,step_angle)


ext=200;
lambdaiext=[lambdai(:,end-ext:end-1)-360 lambdai(:,:) lambdai(:,2:ext)+360];
riext=[ri(:,end-ext:end-1) ri(:,:) ri(:,2:ext)];



delta=min_angle/60:step_angle/60:max_angle/60;

res=zeros(size(delta));


clear angle_min;

progressbar(0);
for j=exclude:(size(ri,1)-exclude)

    for i=1:numel(delta)    
        riis=interp1(lambdaiext(j,:),riext(j,:),...
        lambdai(j,:)+delta(i),'linear');
        res(i)=sum((riis-riref(j,:)).^2);            
    end
    
%     ind=find(res==min(res),1);    
%     angle_min(j-exclude+1)=delta(ind);
%     plot(delta,res); hold on;
    
    p=polyfit(delta,res,2);
    [angle_min(j-exclude+1),fval(j-exclude+1)]=fminsearch(@(x) polyval(p,x),0);
    
%     plot(angle_min(j-exclude+1),fval(j-exclude+1),'x','MarkerSize',4);    
%     unplot(2);    

     progressbar(j/(numel(exclude:(size(ri,1)-exclude))));

end
progressbar(1);


rinew=ri;
progressbar(0)
for j=exclude:(size(ri,1)-exclude)
        rinew(j,:)=interp1(lambdaiext(j,:),riext(j,:),...
        lambdai(j,:)+angle_min(j-exclude+1),'linear');
    progressbar(j/(size(ri,1)-2));
end
progressbar(1);