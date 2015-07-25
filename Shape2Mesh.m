function Shape2Mesh(filename_shape,filename_mesh,nrow)

x = linspace(-1,1,nrow);
y = x;
z = x;

nelems_row = numel(x)-1;
nelems     = nelems_row^3;
[x,y,z]    = ndgrid(x,y,z);
nnodes     = numel(x);

cell_mat   = 0;
cell_type  = 'hex';

s = size(x);

xs = x .* sqrt(1.0 - (y.*y/2.0) - (z.*z/2.0) + (y.*y.*z.*z/3.0));
ys = y .* sqrt(1.0 - (z.*z/2.0) - (x.*x/2.0) + (z.*z.*x.*x/3.0));
zs = z .* sqrt(1.0 - (x.*x/2.0) - (y.*y/2.0) + (x.*x.*y.*y/3.0));

clear x y z

[lon,lat,rs] = cart2sph(xs,ys,zs);

%% Renormalize

TOL      =  1.d-12;

% Open the DSK file for read access.
% We use the DAS-level interface for
% this function.
%
handle = cspice_dasopr( filename_shape );

%
% Begin a forward search through the
% kernel, treating the file as a DLA.
% In this example, it's a very short
% search.
%
[dladsc, found] = cspice_dlabfs( handle );

if ~found
    
    %
    % We arrive here only if the kernel
    % contains no segments. This is
    % unexpected, but we're prepared for it.
    %
    fprintf( 'No segments found in DSK file %s\n', dsk )
    return
    
end

lat = lat(:);
lon = lon(:);

sph = zeros(2,numel(lat));

sph(1,:) = lon;
sph(2,:) = lat;

%
[spoints, ~] = cspice_llgrid_pl02( handle, dladsc, sph );

xsurf = reshape(spoints(1,:),s);
ysurf = reshape(spoints(2,:),s);
zsurf = reshape(spoints(3,:),s);

rsurf = sqrt(xsurf.*xsurf + ysurf.*ysurf + zsurf.*zsurf);

xs = xs.*rsurf;
ys = ys.*rsurf;
zs = zs.*rsurf;

%% Write in file
in = fopen(filename_mesh,'w');

% header

fprintf(in,'%d %d %d %d %d\n',nnodes,nelems,0,0,0);

progressbar(0);

% locations
n = 1:numel(xs);

fprintf(in,'%d %23.16E %23.16E %23.16E\n',[n' xs(:) ys(:) zs(:)]');

% connectivity
ver=zeros(1,8);

cell_num = 1;

for i=1:nelems_row
    for j=1:nelems_row
        for k=1:nelems_row
    
            ver(1) = sub2ind(s, i, j, k);
            ver(2) = sub2ind(s, i+1, j, k);
            ver(3) = sub2ind(s, i+1, j, k+1);
            ver(4) = sub2ind(s, i, j, k+1);
            ver(5) = sub2ind(s, i, j+1, k);
            ver(6) = sub2ind(s, i+1, j+1, k);
            ver(7) = sub2ind(s, i+1, j+1, k+1);
            ver(8) = sub2ind(s, i, j+1, k+1);
    
            fprintf(in,'%d %d %s %d %d %d %d %d %d %d %d\n',...
                cell_num, cell_mat, cell_type,ver);  
            
            cell_num = cell_num + 1;
        end
    end   
    progressbar(i/nelems_row);
end

progressbar(1);

fclose(in);
