function meshStruct = GenerateQuadLayerMesh(lmcosi,lmcosi2,layer_mat,nsq,nl)

s = 100000;

[meshStruct0]=hexMeshSquare([s s],[nsq nsq]);

V = meshStruct0.V;
E = meshStruct0.F;

Eub = find(V(:,2) == s);
Erb = find(V(:,1) == s);

Eb = [Eub; flipud(Erb(1:end-1))];
Vb = V(Eb,:);

[~,latb,~] = cart2sph(Vb(:,1),0,Vb(:,2));

rob = plm2xyz(lmcosi,latb*180/pi,zeros(size(latb)));
[xob,~,zob] = sph2cart(0,latb,rob);

rob2 = plm2xyz(lmcosi2,latb*180/pi,zeros(size(latb)));
[xob2,~,zob2] = sph2cart(0,latb,rob2);

xi = linspacen(Vb(:,1),xob,nl(1));
zi = linspacen(Vb(:,2),zob,nl(1));

[Ei,Vi,~] = surf2patch(xi,zi,zeros(size(xi)),zeros(size(zi))); %Convert to patch data (quadrilateral faces)

Eb = size(Vi,1)-2*nsq:size(Vi,1);

xi2 = linspacen(Vi(Eb,1),xob2,nl(2));
zi2 = linspacen(Vi(Eb,2),zob2,nl(2));

[Ei2,Vi2,~] = surf2patch(xi2,zi2,zeros(size(xi2)),zeros(size(zi2))); %Convert to patch data (quadrilateral faces)

max_ind = max(E(:));

Vtot = [V; Vi];
Etot = [E; Ei+max_ind];
cell_mat = layer_mat(1)+zeros(size(Etot,1),1);

max_ind = max(Etot(:));

Vtot = [Vtot; Vi2];
Etot = [Etot; Ei2+max_ind];
cell_mat = [cell_mat; layer_mat(2)+zeros(size(Ei2,1),1)];

[~,ind1,ind2]=unique(pround(Vtot,5),'rows');
Vtot=Vtot(ind1,:);
Etot=ind2(Etot);

for j=1:size(Etot,1)
    
    x1 = Vtot(Etot(j,1),1);
    y1 = Vtot(Etot(j,1),2);
    
    x2 = Vtot(Etot(j,2),1);
    y2 = Vtot(Etot(j,2),2);
    
    x3 = Vtot(Etot(j,3),1);
    y3 = Vtot(Etot(j,3),2);
    
    x4 = Vtot(Etot(j,4),1);
    y4 = Vtot(Etot(j,4),2);
    
    area = ((x1*y2-y1*x2)+(x2*y3-y2*x3)+(x3*y4-y3*x4)+(x4*y1-y4*x1))/2;
    
    if (area < 0)     
        Etot(j,:) = fliplr(Etot(j,:));       
    end
end

meshStruct.E = Etot;
meshStruct.V = Vtot;
meshStruct.cell_mat = cell_mat;







