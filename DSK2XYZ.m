function [x_surf,y_surf,z_surf] = DSK2XYZ(filename_shape,lat,lon,c)

s = size(lat);

%% Renormalize

TOL      =  1.d-12;

% Open the DSK file for read access.
% We use the DAS-level interface for
% this function.

handle = cspice_dasopr( filename_shape );


% Begin a forward search through the
% kernel, treating the file as a DLA.
% In this example, it's a very short
% search.

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

[spoints, ~] = cspice_llgrid_pl02( handle, dladsc, sph );

x_surf = reshape(spoints(1,:),s)*c;
y_surf = reshape(spoints(2,:),s)*c;
z_surf = reshape(spoints(3,:),s)*c;