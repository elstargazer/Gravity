ccc

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';
[~,shapename,~] = fileparts(shape_filename) ;
shape_full_filename = [shape_folder shape_filename];

filename_mesh = '../CeresFE/FE/mesh_gibbon.inp';

r_sphe  = 470000;
r_core  = 390000;
n_core  = 30;
n_sphe  = 30;

% Shape2Mesh_gibbon(r_sphe, r_core, n_core, n_sphe, ...
%     shape_full_filename, filename_mesh);

RandomShape2Mesh_gibbon(r_sphe, r_core, n_core, n_sphe, ...
    filename_mesh)

L = 60;
filename = '../CeresFE/FE/outer_points.txt';
[sdl,l] = file2spec(filename,L);

% figure; hold on;
% set(gca,'XScale','log');
% set(gca,'YScale','log');

plot(l,sdl,'r-');

xlabel('Degree','FontSize',20);
ylabel('Power','FontSize',20);





%% Generate Jacobi ellipsoid
% T   = 4.6;
% rho = 2169;
% [fh,fval] = HydrostaticStateExact3Ax(r_sphe,T,rho,0.05,0.025);
% [a,b,c] = fr2abc(r_sphe,fh(1),fh(2));
% Ellipsoid2Mesh_gibbon([a b c], r_core, n_core, n_sphe, filename_mesh)

% nrow = 10;
% Shape2Mesh_cubesphere(shape_full_filename,filename_mesh,nrow);
 