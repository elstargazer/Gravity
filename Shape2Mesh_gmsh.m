function Shape2Mesh_gmsh(x,y,z,ver,filename)

lc = 1;
in = fopen(filename,'w');

% Write points
for i=1:numel(x)
    WritePoint_gmsh(in,i,x(i),y(i),z(i),lc)
end

% Write lines and line loops
id_line = 1;
id_lineloop = 1;

for i=1:size(ver,1)
    
    WriteLine_gmsh(in,id_line  ,ver(i,1),ver(i,2));
    WriteLine_gmsh(in,id_line+1,ver(i,2),ver(i,3));
    WriteLine_gmsh(in,id_line+2,ver(i,3),ver(i,4));
    WriteLine_gmsh(in,id_line+3,ver(i,4),ver(i,5));
    
    WriteLineLoop_gmsh(in,id_lineloop,id_line:id_line+3)
    
    id_line = id_line + 4;
    id_lineloop = id_lineloop + 1;
 
end

id_surface = 1;

% Write surface
WriteSurface_gmsh(in,id_surface,1:size(ver,1));

% Write surface loop
WriteSurfaceLoop_gmsh(in,1,1);

% Write volume
WriteVolume_gmsh(in,1,1);

fclose(in);