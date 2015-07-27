% ccc

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';
[~,shapename,~] = fileparts(shape_filename) ;
shape_full_filename = [shape_folder shape_filename];

filename_mesh = '../CeresFE/FE/mesh.inp';
nrow = 10;
Shape2Mesh_cubesphere(shape_full_filename,filename_mesh,nrow);
 
% NTess = 0;
% 
% TR=IcosahedronMesh;
% TR_2=SubdivideSphericalMesh(TR,NTess); 
% 
% FV=TR_2.Triangulation;
% x=TR_2.X(:,1);
% y=TR_2.X(:,2);
% z=TR_2.X(:,3);
% 
% figure; hold on;
% axis equal;
% axis off;
% 
% trisurf(TR_2);
% for i=1:numel(x)
%    
%     plot3(x(i),y(i),z(i),'o','MarkerSize',20);
%     text(1.1*x(i),1.1*y(i),1.1*z(i),num2str(i),'FontSize',20);
%     
% end
% 
% filename_mesh = '../CeresFE/FE/mesh.geo';
% Shape2Mesh_gmsh(x,y,z,FV,filename_mesh);
 
% 