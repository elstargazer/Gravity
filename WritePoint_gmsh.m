function WritePoint_gmsh(in,id,x,y,z,lc)

fprintf(in,'Point(%d) = {%f, %f, %f, %f};\n',id,x,y,z,lc);