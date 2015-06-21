FileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/Corrected/SHAPE_NEW.TXT';

data=load(FileName);

[lambda,fi,r]=cart2sph(data(:,1),data(:,2),data(:,3));

LambdaMin=-180;
LambdaMax=180;
FiMin=-90;
FiMax=90;

LambdaStep=3.6;
FiStep=3.6;


N1=180/FiStep;
N2=360/LambdaStep;

R=zeros(N1,N2);

lambda=lambda*180/pi;
fi=fi*180/pi;



parfor i=1:N1
    for j=1:N2
        
        
        
        Fi1=-90+(i-1)*FiStep;
        Fi2=Fi1+FiStep;         
        Lambda1=-180+(j-1)*LambdaStep;
        Lambda2=Lambda1+LambdaStep;

        
        ConditionLambda=((lambda<Lambda2) & (lambda>Lambda1));
        ConditionFi=((fi<Fi2) & (fi>Fi1));
        
        

        R(i,j)=mean(r(ConditionLambda & ConditionFi));
        
        
    end
end



R=R*1000;

in=fopen('Grid shape Vesta mean2','w');

format='   ';

format_number='%16.9f';

for i=1:N2-1
    format=[format format_number '  '];
end

format=[format format_number '\n'];






for i=1:N1
    
    fprintf(in,format,R(i,:));
    
    
end

fclose(in);