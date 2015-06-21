function MakeLayeredModelMovie(lmcosi1,lmcosi2,lmcosi3)

vidObj = VideoWriter('LayeredModel.avi');
open(vidObj)

FiStep=1;
LambdaStep=1;

Fi=(-90:FiStep:90);

figure('Position',[1 1 1000 1000],'Color',[0 0 0]);
axes('Color',[0 0 0]);

axis([-300 300 -300 300]);

MantleColor=[0.2 0.1 0.6];
ShapeColor=[0.3 0.2 0.9];
CoreColor=[0.05 0.05 0.5];


hold on;
    
%     Lambda=45;

    Lambdai=Fi.*0+Lambda;
    
    WriteXYZ(Lambdai,Fi,Lambdai*0+1,'Crossection1.xyz');
    WriteXYZ(Lambdai+180,Fi,Lambdai*0+1,'Crossection2.xyz');
    
    [ri_shape1_e,~]=plm2xyz(lmcosi1,Fi,Lambdai);
    [ri_shape2_e,~]=plm2xyz(lmcosi2,Fi,Lambdai);
    [ri_shape3_e,~]=plm2xyz(lmcosi3,Fi,Lambdai);
    
    [ri_shape1_w,~]=plm2xyz(lmcosi1,Fi,Lambdai+180);
    [ri_shape2_w,~]=plm2xyz(lmcosi2,Fi,Lambdai+180);
    [ri_shape3_w,~]=plm2xyz(lmcosi3,Fi,Lambdai+180);
    
    [x1_w,y1_w,z1_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape1_w'/1000);
    [x2_w,y2_w,z2_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape2_w'/1000);
    [x3_w,y3_w,z3_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape3_w'/1000);
    
    [x1_e,y1_e,z1_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape1_e'/1000);
    [x2_e,y2_e,z2_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape2_e'/1000);
    [x3_e,y3_e,z3_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape3_e'/1000);
    
    
    x1_w_new=sqrt(x1_w.*x1_w+y1_w.*y1_w);
    x2_w_new=sqrt(x2_w.*x2_w+y2_w.*y2_w);
    x3_w_new=sqrt(x3_w.*x3_w+y3_w.*y3_w);
    
    x1_e_new=sqrt(x1_e.*x1_e+y1_e.*y1_e);
    x2_e_new=sqrt(x2_e.*x2_e+y2_e.*y2_e);
    x3_e_new=sqrt(x3_e.*x3_e+y3_e.*y3_e);
    
    h3=fill(-x1_e_new,z1_e,ShapeColor,'EdgeColor','none');
    h1=fill(x1_w_new,z1_w,ShapeColor,'EdgeColor','none');
    
    h4=fill(-x2_e_new,z2_e,MantleColor,'EdgeColor','none');   
    h2=fill(x2_w_new,z2_w,MantleColor,'EdgeColor','none');
    
    h5=fill(-x3_e_new,z3_e,CoreColor,'EdgeColor','none');   
    h6=fill(x3_w_new,z3_w,CoreColor,'EdgeColor','none');
    
plotplm(lmcosi1,[],[],2,1,[],[],[]);



    currFrame = getframe;
    writeVideo(vidObj,currFrame);

for Lambda=LambdaStep:LambdaStep:180
    
    Lambdai=Fi.*0+Lambda;
    
    [ri_shape1_e,~]=plm2xyz(lmcosi1,Fi,Lambdai);
    [ri_shape2_e,~]=plm2xyz(lmcosi2,Fi,Lambdai);
    [ri_shape3_e,~]=plm2xyz(lmcosi3,Fi,Lambdai);
    
    [ri_shape1_w,~]=plm2xyz(lmcosi1,Fi,Lambdai+180);
    [ri_shape2_w,~]=plm2xyz(lmcosi2,Fi,Lambdai+180);
    [ri_shape3_w,~]=plm2xyz(lmcosi3,Fi,Lambdai+180);
    
    [x1_w,y1_w,z1_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape1_w'/1000);
    [x2_w,y2_w,z2_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape2_w'/1000);
    [x3_w,y3_w,z3_w]=sph2cart(Lambdai/180*pi,Fi/180*pi,ri_shape3_w'/1000);
    
    [x1_e,y1_e,z1_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape1_e'/1000);
    [x2_e,y2_e,z2_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape2_e'/1000);
    [x3_e,y3_e,z3_e]=sph2cart((Lambdai+180)/180*pi,Fi/180*pi,ri_shape3_e'/1000);
    
    
    x1_w_new=sqrt(x1_w.*x1_w+y1_w.*y1_w);
    x2_w_new=sqrt(x2_w.*x2_w+y2_w.*y2_w);
    x3_w_new=sqrt(x3_w.*x3_w+y3_w.*y3_w);
    
    x1_e_new=sqrt(x1_e.*x1_e+y1_e.*y1_e);
    x2_e_new=sqrt(x2_e.*x2_e+y2_e.*y2_e);
    x3_e_new=sqrt(x3_e.*x3_e+y3_e.*y3_e);
    
    
    set(h3,'XData',-x1_e_new,'YData',z1_e);
    set(h1,'XData',x1_w_new,'YData',z1_w);
    set(h4,'XData',-x2_e_new,'YData',z2_e);
    set(h2,'XData',x2_w_new,'YData',z2_w);
    set(h5,'XData',-x3_e_new,'YData',z3_e);
    set(h6,'XData',x3_w_new,'YData',z3_w);
   
    text(-0,-0,['Lambda=' num2str(Lambda)]);

    
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
    
    unplot
    
    
end


close(vidObj);