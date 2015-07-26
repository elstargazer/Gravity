function WriteLine_gmsh(in,id,p1,p2)

fprintf(in,'Line(%d) = {%d, %d};\n',id,p1,p2);