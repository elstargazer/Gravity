
figure('Color','k','Position',[1 1 1000 1000]); hold on;
%set(gca,'FontSize',20,'Position', [0 0 1 1]);

[lambda_grid,fi_grid,r_grid]=cart2sph(x_grid,y_grid,z_grid);

h=surf(x_grid*1000,y_grid*1000,z_grid*1000,H_grid);
set(h,'BackFaceLighting','lit');

colormap parula;
material([0 1 1]);
lighting phong
shading interp
axis equal tight off
box off

view(0,0);

lambda_light=10+180;
fi_light=0;

[x_light,y_light,z_light]=sph2cart(lambda_light/180*pi,...
    fi_light/180*pi,1);

light_handle=light('Style','infinite',...
'Position',[x_light y_light z_light]);

MovieFileName='~/Dawn/Figures/Ceres/CeresRotation4.avi';
vidObj = VideoWriter(MovieFileName);
open(vidObj);
axis vis3d

factor=0:0.025:3.5;

for i=1:numel(factor)
    
    r_grid_mod=r_grid+factor(i)*H_grid;
    
    [x_grid_mod,y_grid_mod,z_grid_mod]=sph2cart(lambda_grid,...
        fi_grid,r_grid_mod);
     
    set(h,'XData',x_grid_mod*1000,...
    'YData',y_grid_mod*1000,...
    'ZData',z_grid_mod*1000)
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end

ang=0:1:720;

for i=1:numel(ang)
    
    lambda_light=ang(i)+10+180;
    fi_light=0;
    
    [x_light,y_light,z_light]=sph2cart(lambda_light/180*pi,...
        fi_light/180*pi,1);
    
    view(ang(i),0);
    set(light_handle,'Position',[x_light y_light z_light]);
    drawnow
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end


factor=fliplr(factor);

for i=1:numel(factor)
    
    r_grid_mod=r_grid+factor(i)*H_grid;
    
    [x_grid_mod,y_grid_mod,z_grid_mod]=sph2cart(lambda_grid,...
        fi_grid,r_grid_mod);
     
    set(h,'XData',x_grid_mod*1000,...
    'YData',y_grid_mod*1000,...
    'ZData',z_grid_mod*1000)
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end


close(vidObj);
