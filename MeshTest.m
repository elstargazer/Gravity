% ccc

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';
[~,shapename,~] = fileparts(shape_filename) ;
shape_full_filename = [shape_folder shape_filename];

filename_mesh = '../CeresFE/FE/mesh.inp';

nrow = 60;

Shape2Mesh(shape_full_filename,filename_mesh,nrow);