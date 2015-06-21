ccc
tic
%InFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POSTLAMO/VEST64_DTM.raw';
%OutFileName='VestaPOSTLAMOshape_geod_elev.raw';

InFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/Vesta2014/VEST64_DTM.raw';
OutFileName='Vesta20140513shape_geoc_elev64.raw';

% POSTLAMO parameters
% par1=250;
% par2=0.3809711379D+02;
% par3=0.1234695905D-02;

%HAMO2 parameters
par1=250.00000;
par2=0.3806461345D+02;
par3=0.1234327396D-02;


[~,~,R]=ReadRawShapeModel(InFileName,par1,par2,par3,[23041 11521]);

1

[offset,scale]=WriteRawShapeModel(OutFileName,R);

offset
scale

NaN_Fraction=sum(sum(isnan(R)))/numel(R)


2
toc