function WriteSurface_gmsh(in,id,pl)

fprintf(in,'Plane Surface(%d) = {',id);

for i=1:numel(pl)-1 
    fprintf(in,'%d, ',pl(i));         
end
fprintf(in,'%d};\n',pl(end));      

