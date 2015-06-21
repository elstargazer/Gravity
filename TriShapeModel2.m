tic

%% Input parameters

NTess=7; % number of tesselations
MaxDegreeTopo=150; % max SH degree of topo

%% Make Triangular Mesh
TR=IcosahedronMesh;

TR_2=SubdivideSphericalMesh(TR,NTess); 
% figure, h=trimesh(TR_2); set(h,'EdgeColor','b'), axis equal

FV=TR_2.Triangulation;
x_t=TR_2.X(:,1);
y_t=TR_2.X(:,2);
z_t=TR_2.X(:,3);

[lambda_t,fi_t,~]=cart2sph(x_t,y_t,z_t);

%% Load SH shape model
% load Vesta20130522shape_sh1800.mat
% lmcosi_shape=TruncateGravityModel(lmcosi_shape,MaxDegreeTopo,1);
%% Load grid shape model

filename='grd/Vesta20140513shape_geoc_elev_3m.grd';
[lambda_s,fi_s,r_s]=ReadGRD(filename);

N=numel(lambda_t);

r_s_mod=r_s;

r_s_mod((1:size(imi,1))+1500,1:size(imi,2))=...
    r_s((1:size(imi,1))+1500,1:size(imi,2))+imi*14;

tic
r=griddata(lambda_s,fi_s,r_s_mod,lambda_t(1:N)*180/pi,fi_t(1:N)*180/pi,'nearest');
toc

%% Shape factor 
% Degree=lmcosi_shape(:,1);
% figure; hold on;
% 
% xlim([3 max(Degree(:))]);
% ylim([0 20]);
% 
% gi=ginput(3);
% 
% f=interp1(gi(:,1),gi(:,2),Degree,'cubic');
% 
% plot(Degree,f,'.');
% 
% lmcosi_shape(:,3)=lmcosi_shape(:,3).*f;
% lmcosi_shape(:,4)=lmcosi_shape(:,4).*f;

%% Computing radius at verteces

% [r,~,~,~]=plm2xyz(lmcosi_shape,fi_t*180/pi,lambda_t*180/pi,[],[]);


%% Computing x y z

[xi,yi,zi]=sph2cart(lambda_t,fi_t,r);


%% Plotting shape model
figure

tri_fig=trisurf(FV,xi,yi,zi,r); shading interp
colorbar; axis equal

toc
save('VestaTriShape4.mat','xi','yi','zi','FV');

