 ccc

filename='/Users/antonermakov/Downloads/EAR_A_I0037_5_BENNUSHAPE_V1_0/data/101955bennu.tab';
[x,y,z,Ver1,Ver2,Ver3]=ReadTriShapeModel2(filename);



figure('color','k','Position',[1 1 2560 1440]);
hold on;
% set(gca,'FontSize',20);

r=sqrt(x.*x+y.*y+z.*z);

trisurf([Ver1' Ver2' Ver3'],x,y,z,r);

lighting phong
shading interp
axis equal tight off
box off
hold on;
light_handle=light('Style','infinite');
set(h,'BackFaceLighting','lit');
colormap jet;
material([0 0.9 0]);



% plot_star=plot3(x_star,y_star,z_star,'.w','markersize',1);
% 
% r_star=[x_star,y_star,z_star];
% r=1500000;
% axis([-r r -r r -r r]);


set(gca, 'Position', [0 0 1 1]);
% T = get(gca,'tightinset');
% set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

view(0,-75)
