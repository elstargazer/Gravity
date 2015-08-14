function Write_ucd(meshStruct,filename_mesh,cell_type)

nnodes = size(meshStruct.V,1);
nelems = size(meshStruct.E,1);

in = fopen(filename_mesh,'w');

fprintf(in,'%d %d %d %d %d\n',nnodes,nelems,0,0,0);

% locations
n = 1:nnodes;
fprintf(in,'%d %13.6E %13.6E %13.6E\n',...
    [n' meshStruct.V(:,1) meshStruct.V(:,2) meshStruct.V(:,3)]');

cell_num = 1;

switch cell_type
    case 'hex'
        ind_str = '%d %d %d %d %d %d %d %d\n';     
    case 'quad'
        ind_str = '%d %d %d %d\n';              
end

for i=1:nelems   
    fprintf(in,['%d %d %s ' ind_str],...
        cell_num, meshStruct.cell_mat(i),...
        cell_type,meshStruct.E(i,:));
    cell_num = cell_num + 1;
end

fclose(in);


