function meshStruct = GenerateQuadLayerMesh(lmcosi,lmcosi2,layer_mat,nsq,nl)

[meshStruct0]=hexMeshSquare([1 1],[nsq nsq]);

V = meshStruct0.V;
E = meshStruct0.F;

Eub = find(V(:,2) == 1);
Erb = find(V(:,1) == 1);

Eb = [Eub; flipud(Erb(1:end))];
Vb = V(Eb,:);

[~,latb,~] = cart2sph(Vb(:,1),0,Vb(:,2));
rob = plm2xyz(lmcosi,latb*180/pi,zeros(size(latb)));
[xob,~,zob] = sph2cart(0,latb,rob);

rob2 = plm2xyz(lmcosi2,latb*180/pi,zeros(size(latb)));
[xob2,~,zob2] = sph2cart(0,latb,rob2);

xi = linspacen(Vb(:,1),xob,nl);
zi = linspacen(Vb(:,2),zob,nl);

[Ei,Vi,~] = surf2patch(xi,zi,zeros(size(xi)),zeros(size(zi))); %Convert to patch data (quadrilateral faces)

Eb = size(Vi,1)-2*nsq-1:size(Vi,1);

xi2 = linspacen(Vi(Eb,1),xob2,nl);
zi2 = linspacen(Vi(Eb,2),zob2,nl);

[Ei2,Vi2,~] = surf2patch(xi2,zi2,zeros(size(xi2)),zeros(size(zi2))); %Convert to patch data (quadrilateral faces)

max_ind = max(E(:));

Vtot = [V; Vi];
Etot = [E; Ei+max_ind];
cell_mat = layer_mat(1)+zeros(size(Etot,1),1);

max_ind = max(Etot(:));

Vtot = [Vtot; Vi2];
Etot = [Etot; Ei2+max_ind];
cell_mat = [cell_mat; layer_mat(2)+zeros(size(Ei2,1),1)];

meshStruct.E = Etot;
meshStruct.V = Vtot;
meshStruct.cell_mat = cell_mat;





