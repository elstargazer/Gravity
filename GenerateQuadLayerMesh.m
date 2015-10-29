function meshStruct = GenerateQuadLayerMesh(lmcosi,nsq,nl)

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


xi = linspacen(Vb(:,1),xob,nl);
zi = linspacen(Vb(:,2),zob,nl);

[Ei,Vi,~] = surf2patch(xi,zi,zeros(size(xi)),zeros(size(zi))); %Convert to patch data (quadrilateral faces)

max_ind = max(E(:));

Vtot = [V; Vi];
Etot = [E; Ei+max_ind];

meshStruct.E = Etot;
meshStruct.V = Vtot;
meshStruct.cell_mat = zeros(size(Etot,1),1);





