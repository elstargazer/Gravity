function Shape2Mesh_cubesphere(filename_shape,filename_mesh,nrow)



%% Parameters 

npts_shell = 10;
cube_len   = 100000;

%% 
xc = linspace(-1,1,nrow)*cube_len;
yc = xc;
zc = xc;

nelems_row = numel(xc)-1;
nelems     = nelems_row^3;
[xc,yc,zc] = ndgrid(xc,yc,zc);
nnodes     = numel(xc);

cell_mat   = 0;
cell_type  = 'hex';

% extract boundary surfaces of the inner cube

xc_bnry_1 = squeeze(xc(1,:,:));
yc_bnry_1 = squeeze(yc(1,:,:));
zc_bnry_1 = squeeze(zc(1,:,:));

xc_bnry_2 = squeeze(xc(end,:,:));
yc_bnry_2 = squeeze(yc(end,:,:));
zc_bnry_2 = squeeze(zc(end,:,:));

xc_bnry_3 = squeeze(xc(:,1,:));
yc_bnry_3 = squeeze(yc(:,1,:));
zc_bnry_3 = squeeze(zc(:,1,:));

xc_bnry_4 = squeeze(xc(:,end,:));
yc_bnry_4 = squeeze(yc(:,end,:));
zc_bnry_4 = squeeze(zc(:,end,:));

xc_bnry_5 = xc(:,:,1);
yc_bnry_5 = yc(:,:,1);
zc_bnry_5 = zc(:,:,1);

xc_bnry_6 = xc(:,:,end);
yc_bnry_6 = yc(:,:,end);
zc_bnry_6 = zc(:,:,end);

% compute lat and lon of the radial lines

[lon_bnry_1,lat_bnry_1,~] = cart2sph(xc_bnry_1,yc_bnry_1,zc_bnry_1);
[lon_bnry_2,lat_bnry_2,~] = cart2sph(xc_bnry_2,yc_bnry_2,zc_bnry_2);
[lon_bnry_3,lat_bnry_3,~] = cart2sph(xc_bnry_3,yc_bnry_3,zc_bnry_3);
[lon_bnry_4,lat_bnry_4,~] = cart2sph(xc_bnry_4,yc_bnry_4,zc_bnry_4);
[lon_bnry_5,lat_bnry_5,~] = cart2sph(xc_bnry_5,yc_bnry_5,zc_bnry_5);
[lon_bnry_6,lat_bnry_6,~] = cart2sph(xc_bnry_6,yc_bnry_6,zc_bnry_6);

% getting coordinates on the outer surface
[x_surf_1,y_surf_1,z_surf_1] = DSK2XYZ(filename_shape,lat_bnry_1,lon_bnry_1,1000);
[x_surf_2,y_surf_2,z_surf_2] = DSK2XYZ(filename_shape,lat_bnry_2,lon_bnry_2,1000);
[x_surf_3,y_surf_3,z_surf_3] = DSK2XYZ(filename_shape,lat_bnry_3,lon_bnry_3,1000);
[x_surf_4,y_surf_4,z_surf_4] = DSK2XYZ(filename_shape,lat_bnry_4,lon_bnry_4,1000);
[x_surf_5,y_surf_5,z_surf_5] = DSK2XYZ(filename_shape,lat_bnry_5,lon_bnry_5,1000);
[x_surf_6,y_surf_6,z_surf_6] = DSK2XYZ(filename_shape,lat_bnry_6,lon_bnry_6,1000);

% interpolating values in the shell
x_shell_1 = linspaceNDim(xc_bnry_1, x_surf_1, npts_shell);
y_shell_1 = linspaceNDim(yc_bnry_1, y_surf_1, npts_shell);
z_shell_1 = linspaceNDim(zc_bnry_1, z_surf_1, npts_shell);

x_shell_2 = linspaceNDim(xc_bnry_2, x_surf_2, npts_shell);
y_shell_2 = linspaceNDim(yc_bnry_2, y_surf_2, npts_shell);
z_shell_2 = linspaceNDim(zc_bnry_2, z_surf_2, npts_shell);

x_shell_3 = linspaceNDim(xc_bnry_3, x_surf_3, npts_shell);
y_shell_3 = linspaceNDim(yc_bnry_3, y_surf_3, npts_shell);
z_shell_3 = linspaceNDim(zc_bnry_3, z_surf_3, npts_shell);

x_shell_4 = linspaceNDim(xc_bnry_4, x_surf_4, npts_shell);
y_shell_4 = linspaceNDim(yc_bnry_4, y_surf_4, npts_shell);
z_shell_4 = linspaceNDim(zc_bnry_4, z_surf_4, npts_shell);

x_shell_5 = linspaceNDim(xc_bnry_5, x_surf_5, npts_shell);
y_shell_5 = linspaceNDim(yc_bnry_5, y_surf_5, npts_shell);
z_shell_5 = linspaceNDim(zc_bnry_5, z_surf_5, npts_shell);

x_shell_6 = linspaceNDim(xc_bnry_6, x_surf_6, npts_shell);
y_shell_6 = linspaceNDim(yc_bnry_6, y_surf_6, npts_shell);
z_shell_6 = linspaceNDim(zc_bnry_6, z_surf_6, npts_shell);


% plot points

ccj = jet(6);

figure; hold on;
axis equal;

plot3(x_shell_1(:),y_shell_1(:),z_shell_1(:),'.','Color',ccj(1,:));
plot3(x_shell_2(:),y_shell_2(:),z_shell_2(:),'.','Color',ccj(2,:));
plot3(x_shell_3(:),y_shell_3(:),z_shell_3(:),'.','Color',ccj(3,:));
plot3(x_shell_4(:),y_shell_4(:),z_shell_4(:),'.','Color',ccj(4,:));
plot3(x_shell_5(:),y_shell_5(:),z_shell_5(:),'.','Color',ccj(5,:));
plot3(x_shell_6(:),y_shell_6(:),z_shell_6(:),'.','Color',ccj(6,:));

plot3(xc(:),yc(:),zc(:),'.','Color','k');

%% Write in file
in = fopen(filename_mesh,'w');

% header
fprintf(in,'%d %d %d %d %d\n',nnodes,nelems,0,0,0);
progressbar(0);

% locations
n = 1:numel(xc);

fprintf(in,'%d %23.16E %23.16E %23.16E\n',[n' xc(:) yc(:) zc(:)]');

% connectivity
ver=zeros(1,8);

cell_num = 1;

for i=1:nelems_row
    for j=1:nelems_row
        for k=1:nelems_row
    
            ver(1) = sub2ind(s, i, j, k);
            ver(2) = sub2ind(s, i+1, j, k);
            ver(3) = sub2ind(s, i+1, j, k+1);
            ver(4) = sub2ind(s, i, j, k+1);
            ver(5) = sub2ind(s, i, j+1, k);
            ver(6) = sub2ind(s, i+1, j+1, k);
            ver(7) = sub2ind(s, i+1, j+1, k+1);
            ver(8) = sub2ind(s, i, j+1, k+1);
    
            fprintf(in,'%d %d %s %d %d %d %d %d %d %d %d\n',...
                cell_num, cell_mat, cell_type,ver);  
            
            cell_num = cell_num + 1;
        end
    end   
    progressbar(i/nelems_row);
end

progressbar(1);

fclose(in);