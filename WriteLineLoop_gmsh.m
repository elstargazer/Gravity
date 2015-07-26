function WriteLineLoop_gmsh(in,id,lns)

fprintf(in,'Line Loop(%d)',id);

for i=1:numel(lns)    
   fprintf(in,'%d, ',lns(i));
end

fprintf(in,'};\n');
