function Shape2Mesh_gmsh(x,y,z,ver,filename)

lc = 0.5; % characteristic length

id_point         = 1;
id_line          = 1;
id_lineloop      = 1;
id_surface       = 1;
id_surface_loop  = 1;
id_volume        = 1;

% open file
in = fopen(filename,'w');

[path,~,~] = fileparts(filename) ;

debugp_filename = [path '/debugp.pos'];
debugr_filename = [path '/debugr.pos'];

in_p = fopen(debugp_filename,'w');
in_r = fopen(debugr_filename,'w');

fclose(in_p);
fclose(in_r);

fprintf(in,'Mesh.RandomFactor=1e-4;\n\n');

%% Write points
for i=1:numel(x)
    WritePoint_gmsh(in,id_point,x(i),y(i),z(i),lc)
    id_point = id_point + 1;
end

%% Write lines and line loops
s = size(ver);

id_line_vec     =  zeros(1,s(2));
id_surface_vec  =  zeros(1,s(1));

for i=1:s(1)
    
    % write lines   
    for j=1:s(2)-1
        WriteLine_gmsh(in,id_line,ver(i,j),ver(i,j+1));
        id_line_vec(j) = id_line;
        id_line = id_line + 1;
    end  
    WriteLine_gmsh(in,id_line,ver(i,j+1),ver(i,j-s(2)+2));
    id_line_vec(j+1) = id_line;
    id_line          = id_line + 1;
    
    % write line loop
    WriteLineLoop_gmsh(in,id_line,id_line_vec) 
    % write surface 
    WriteSurface_gmsh(in,id_surface,id_line);
    id_surface_vec(i) = id_surface;
    id_surface        = id_surface + 1;
    id_line = id_line + 1;
end

%% Write surface loop
WriteSurfaceLoop_gmsh(in,id_surface,id_surface_vec);

%% Write volume
WriteVolume_gmsh(in,id_volume,id_surface);

fclose(in);