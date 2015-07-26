function WriteSurface_gmsh(in,id,pl)

fprintf('Plane Surface(%d) = {',id);

for i=1:numel(pl)   
    fprintf('%d',pl(i));      
end

fprintf(in,'};\n');