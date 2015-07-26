function WriteVolume_gmsh(in,id,sfs_loops)

fprintf('Volume (%d) = {',id);

for i=1:numel(pl)   
    fprintf('%d',sfs_loops(i));      
end

fprintf(in,'};\n');