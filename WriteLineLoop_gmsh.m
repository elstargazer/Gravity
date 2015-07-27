function WriteLineLoop_gmsh(in,id,lns)

fprintf(in,'Line Loop(%d) = {',id);

for i=1:numel(lns)-1  
   fprintf(in,'%d, ',lns(i));
end

fprintf(in,'%d};\n',lns(end));


