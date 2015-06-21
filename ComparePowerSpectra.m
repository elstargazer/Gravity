


% List of shapes
ShapeListFile = 'ShapeList.txt';

% Max spherical harmonics degree

L = 20;
step = 1;

in = fopen(ShapeListFile,'r');
ShapeFile='shape';

% Plot
fig = figure('Position',[1 1 1300 1000]);
set(gca,'FontSize',20); hold on;
set(gca,'YScale','log');

while (ShapeFile~=-1)
    
    ShapeFile = fgetl(in);
    [x,y,z,Ver1,Ver2,Ver3]=ReadTriShapeModelGaskell(ShapeFile);   
    [lambdai,fii,ri]=Tri2Grid([x; y; z],[Ver1; Ver2; Ver3],step);
%     lmcosi_shape = xyz2plm(ri,L);
%     [sdl,l]=plm2spec(lmcosi_shape);
%     plot(l,sdl,'-k','LineWidth',3);
    
end

fclose(in);