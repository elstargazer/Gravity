% ccc


%% Getting triangular mesh
% TriShapeModel2

%% Scaling
Scale=1/10e6*1000;

ver=1000*Scale*[xi yi zi];
 
% cond = ver(:,3)>-39;
% ver_to_delete = find(cond);
% 
% del = zeros(1,size(FV,1));
% 
% progressbar(0)
% parfor i=1:size(FV,1)   
%     del(i) = ismember(FV(i,1),ver_to_delete) ||...
%         ismember(FV(i,2),ver_to_delete) || ...
%         ismember(FV(i,3),ver_to_delete); 
% end
% progressbar(1);
% 
% FV2=FV;
% ver2=ver;
% 
% progressbar(0);
% j=1;
% for i=size(FV2,1):-1:1
%     
%     if (del(i)==1)
%         FV2(i,:)=[];
%     end
%     j=j+1;
%     progressbar(j/size(FV2,1));
% end
% progressbar(1);


% j=1;
% progressbar(0);
% ni=size(ver2,1);
% for i=ni:-1:1
%     if (cond(i)==1)
%         ver2(i,:)=[];
%     end
%     j=j+1;
%     progressbar(j/ni);
% end
% progressbar(1);

A.faces=FV;
A.vertices=ver;

figure; hold on; axis equal;
tri_fig=trisurf(FV,ver(:,1), ver(:,2), ver(:,3),r); shading interp

xrange=range(ver(:,1))
yrange=range(ver(:,2))
zrange=range(ver(:,3))

[xrange yrange zrange]

% filename_ascii=['VestaShapeModel_' num2str(NTess) '_' ...
%     num2str(MaxDegreeTopo) '_ascii.stl'];

filename_ascii=['VestaShapeModel_Dinara2_' num2str(NTess) '_' ...
     '3m' '_ascii.stl'];
 
 filename_bin=['VestaShapeModel_Dinara2_' num2str(NTess) '_' ...
     '3m' '_bin.stl'];

% filename_bin='VestaShapeModel_3_20_bin.stl';

stlwrite(filename_ascii, A,'MODE','ascii');
stlwrite(filename_bin, A,'MODE','binary');

%% STL from gridded file

% filename_ascii='Vesta20130522shape_geoc_elev_12m_gridline.stl';
% 
% 
% filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
% [lambda_s,fi_s,r_s]=ReadGRD(filename);
% 
% [x,y,z]=sph2cart(lambda_s/180*pi,fi_s/180*pi,r_s*1000*Scale);

% stlwrite(filename_ascii, x,y,z,'MODE','ascii','TRIANGULATION','x');