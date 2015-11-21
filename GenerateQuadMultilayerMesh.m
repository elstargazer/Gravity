function meshStruct = GenerateQuadMultilayerMesh(lmcosi_cell,layer_mat,nsq,nl)

n_layer = numel(lmcosi_cell);

layer_mat = flipud(layer_mat);

%% Create inner cube mesh
s = 50000;

[meshStruct0]=hexMeshSquare([s s],[nsq nsq]);

V = meshStruct0.V;
E = meshStruct0.F;

Eub = find(V(:,2) == s);
Erb = find(V(:,1) == s);

Eb = [Eub; flipud(Erb(1:end-1))];
Vb = V(Eb,:);

% latitude of the box boundary
[~,latb,~] = cart2sph(Vb(:,1),0,Vb(:,2));

% inner layer coordinates
xil = Vb(:,1);
zil = Vb(:,2);

Etot = E;
Vtot = V;
cell_mat_tot = layer_mat(1) + zeros(size(Etot,1),1);

for i=1:n_layer
    
    % outer surface coordinates
    rol = plm2xyz(lmcosi_cell{i},latb*180/pi,zeros(size(latb)));
    [xol,~,zol] = sph2cart(0,latb,rol);
    
    % find points between layers uniformly spaced
    xi = linspacen(xil,xol,nl(i));
    zi = linspacen(zil,zol,nl(i));
    
    % convert uniformly spaced points to patches (quads)
    [Ei,Vi,~] = surf2patch(xi,zi,zeros(size(xi)),zeros(size(zi))); %Convert to patch data (quadrilateral faces)
    
    % indices of the boundary layer
    Eb = size(Vi,1)-2*nsq:size(Vi,1);
    
    max_ind = max(Etot(:));
    
    Vtot = [Vtot; Vi];
    Etot = [Etot; Ei+max_ind];
    cell_mat_tot = [cell_mat_tot; layer_mat(i)+zeros(size(Ei,1),1)];
    
    xil = Vi(Eb,1);
    zil = Vi(Eb,2);
     
end

% check if there any duplicate vertices
[~,ind1,ind2]=unique(pround(Vtot,5),'rows');
Vtot=Vtot(ind1,:);
Etot=ind2(Etot);

%% invert if volume is negative

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

% create a mesh structure
meshStruct.E = Etot;
meshStruct.V = Vtot;
meshStruct.cell_mat = cell_mat_tot;







