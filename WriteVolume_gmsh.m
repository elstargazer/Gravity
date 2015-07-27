function WriteVolume_gmsh(in,id,sfs_loops)

fprintf(in,'Volume (%d) = {',id);

for i=1:numel(sfs_loops)-1   
    fprintf(in,'%d, ',sfs_loops(i));  
end

fprintf(in,'%d};\n',sfs_loops(end));  
