ccc
tic
%InFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POSTLAMO/VEST64_DTM.raw';
%OutFileName='VestaPOSTLAMOshape_geod_elev.raw';

InFileName='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POST_VESTA/VEST64_DTM.raw';
OutFileName='VestaPOSTVESTAshape_geod_elev.raw';

% POSTLAMO parameters
% par1=250;
% par2=0.3809711379D+02;
% par3=0.1234695905D-02;

%HAMO2 parameters
par1=250;
par2=0.3806459373D+02;
par3=0.1234318295D-02;

[lambda,fi,R]=ReadRawShapeModel(InFileName,par1,par2,par3);

1

lambda=lambda/180*pi;
fi=fi/180*pi;
[x,y,z]=sph2cart(lambda,fi,R);

2

clear fi lambda R
% a=280.857; %% Topography ellipsoid
% c=226.421;

a=280.871;
c=226.505; %%  ellipsoid

e=Eccentricity(a,c);
[B,L,H]=XYZ2BLH(x,y,z,a,e);

3

clear x y z

[offset,scale]=WriteRawShapeModel(OutFileName,H);

NaN_Fraction=sum(sum(isnan(H)))/numel(H)


4
toc