function Shape2Mesh_gibbon(r_sph, r_core, ncore, nmantle, shape_file, mesh_file)


%% Creating a solid hexahedral mesh  

%Control settings
cPar.sphereRadius=r_sph;
cPar.coreRadius=r_core;
cPar.numElementsMantel=nmantle; 
cPar.numElementsCore=ncore; 
cPar.makeHollow=0;

%Creating mesh
tic
% [meshStruct]=hexMeshSphere(cPar);
[meshStruct]=hexMeshShape(cPar,shape_file);

toc

%Access ouput
E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

%% Write to ucd

Write_ucd(meshStruct,mesh_file)

%% Plotting mesh 

% Plot settings
figColor     = 'w'; 
figColorDef  = 'white';
fontSize     = 15;
faceAlpha1   = 1;
edgeColor    = 0.25*ones(1,3);
edgeWidth    = 2;

% plotting

hf=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('The hexahedral mesh sphere','FontSize',fontSize);
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
L=YE>mean(Y);
[Fs,~]=element2patch(E(L,:),[],'hex8');
patch('Faces',Fs,'Vertices',V,'FaceColor','b',...
    'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth/2,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

drawnow; 