function Write_ucd(meshStruct,filename_mesh)

cell_mat   = 0;
cell_type  = 'hex';

E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

nnodes = size(V,1);
nelems = size(E,1);

in = fopen(filename_mesh,'w');

fprintf(in,'%d %d %d %d %d\n',nnodes,nelems,0,0,0);

% locations
n = 1:nnodes;
fprintf(in,'%d %13.6E %13.6E %13.6E\n',...
    [n' V(:,1) V(:,2) V(:,3)]');

cell_num = 1;

for i=1:nelems   
    fprintf(in,'%d %d %s %d %d %d %d %d %d %d %d\n',...
        cell_num, cell_mat, cell_type,E(i,:));   
end



fclose(in);


