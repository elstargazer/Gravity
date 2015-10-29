ccc

filename_mesh = 'quad_new_mesh.inp';
nsq = 40;
nl  = 10;


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
      
meshStruct = GenerateQuadLayerMesh(lmcosi,nsq,nl);

V = meshStruct.V;
E = meshStruct.E;

figure; hold on;
axis equal;

for j=1:size(E,1)
    
    line([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
    
end

 Write_ucd(meshStruct,filename_mesh,'quad');




