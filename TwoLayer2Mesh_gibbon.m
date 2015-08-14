function TwoLayer2Mesh_gibbon(lmcosi1, lmcosi2, n_1, n_2, n_cube, r_cube, mesh_file)


%% Creating a solid hexahedral mesh  

r1 = lmcosi1(1,3);
r2 = lmcosi2(1,3);

%Control settings
cPar.cubeRadius        = r_cube;
cPar.numElementsCube   = n_cube; 
cPar.numElementsMantel = n_2; 
cPar.numElementsCore   = n_1; 

%Creating sphere
tic
% [meshStruct]=hexMeshSphere(cPar);
[meshStruct]=hexTwoLayerShape(cPar,lmcosi1,lmcosi2);
toc

%Access ouput
E  = meshStruct.E; %The elements 
V  = meshStruct.V; %The vertices
Fb = meshStruct.Fb; %The boundary faces
cell_mat = meshStruct.cell_mat;


meshStruct

%% Write to ucd

Write_ucd(meshStruct,mesh_file,'hex');

%% Plotting mesh 

% Plot settings
figColor    = 'w'; 
figColorDef = 'white';
fontSize    = 15;
faceAlpha1  = 1;
edgeColor   = 0.25*ones(1,3);
edgeWidth   = 2;

% plotting

hf=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('The hexahedral mesh of a random shape','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','r','FaceAlpha',faceAlpha1,...
    'lineWidth',edgeWidth/2,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

subplot(1,2,2);
title('Cut-view of the mesh','FontSize',fontSize);

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>0;
cell_mat = cell_mat(L);
[Fs,~]=element2patch(E(L,:),[],'hex8');

patch('Faces',Fs,'Vertices',V,'FaceColor','b',...
    'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth/2,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

drawnow; 



