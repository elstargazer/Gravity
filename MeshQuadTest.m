ccc

filename_mesh = 'quad_new_mesh.inp';
nsq = 20;
nl  = 20;


lmcosi = [0 0 3 0;
    1 0 0 0;
    1 1 0 0;
    2 0 0 0;
    2 1 0 0;
    2 2 0 0;
    3 0 0.2 0;
    3 1 0 0;
    3 2 0 0;
    3 3 0 0];

lmcosi2 = [0 0 6 0;
    1 0 0 0;
    1 1 0 0;
    2 0 0 0;
    2 1 0 0;
    2 2 0 0;
    3 0 -0.4 0;
    3 1 0 0;
    3 2 0 0;
    3 3 0 0];

layer_mat = [0 1];

meshStruct = GenerateQuadLayerMesh(lmcosi,lmcosi2,layer_mat,nsq,nl);

V = meshStruct.V;
E = meshStruct.E;
cell_mat = meshStruct.cell_mat;

figure; hold on;
axis equal;

for j=1:size(E,1)
    
    if cell_mat(j) == 0;
        Color = [1 0 0];
    elseif cell_mat(j) == 1;
        Color = [0 0 1];
    else
        
    end
    
    patch([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],...
        Color);
    
    line([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
end

Write_ucd(meshStruct,filename_mesh,'quad');




