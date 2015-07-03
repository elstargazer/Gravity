function [x,y,z]=ReadSPC(filename,A,opt)

TOL      =  1.d-12;

% Open the DSK file for read access.
% We use the DAS-level interface for
% this function.
%
handle = cspice_dasopr( filename );

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

%
% If we made it this far, DLADSC is the
% DLA descriptor of the first segment.

if (opt == 'grid')
    
    step = A;
        
    lat=(90:-step:-90) * cspice_rpd();
    lon=(0:step:360) * cspice_rpd();
    
    [lati, loni] = meshgrid(lat, lon);
    
    s = size( lati );
    
    lati = lati(:);
    loni = loni(:);
    
    grid = zeros(2,numel(lati));
    
    grid(1,:) = loni;
    grid(2,:) = lati;
    
    %
    [spoints, ~] = cspice_llgrid_pl02( handle, dladsc, grid );
    
    x = spoints(1,:);
    y = spoints(2,:);
    z = spoints(3,:);
    
    x = reshape(x,s);
    y = reshape(y,s);
    z = reshape(z,s);
    
elseif (opt == 'rand')
    
    Npoints = A;
    
    [fii,lambdai]=GenerateRandomSphericalCoord(Npoints);
   
    [spoints, ~] = cspice_llgrid_pl02( handle, dladsc, [lambdai'; fii'] );
    
    x = spoints(1,:);
    y = spoints(2,:);
    z = spoints(3,:);
    
end
    
cspice_dascls( handle );

    
    
    
    
    
    