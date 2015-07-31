function meshStruct = Read_ucd(filename)

in = fopen(filename);

% Read header
C = textscan(in,'%d %d %d %d %d\n',1);
nnodes = (C{1});
nelems = (C{2});

% Read body
data = textscan(in,'%d %f %f %f\n',nnodes);
meshStruct.V = ...
    [(data{2}) (data{3}) (data{4})];
data = textscan(in,'%d %d %s %d %d %d %d\n',nelems);
meshStruct.E = [(data{4}) (data{5}) ...
    (data{6}) (data{7})];

fclose(in);