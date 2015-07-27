function WriteSurfaceLoop_gmsh(in,id,sfs)

fprintf(in,'Surface Loop(%d) = {',id);

for i=1:numel(sfs)-1    
   fprintf(in,'%d, ',sfs(i));
end

fprintf(in,'%d};\n',sfs(end));