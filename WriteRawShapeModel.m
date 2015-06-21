function [offset,scale]=WriteRawShapeModel(FileName,H)

in=fopen(FileName,'w');


H_min=min(min(H));
H_max=max(max(H));

H_range=H_max-H_min;

scale=H_range/65535;
offset=H_min;

H=(H-offset)/scale;

fwrite(in,H,'uint16','b');

fclose(in);