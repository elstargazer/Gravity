function [lambdai,fii,ri]=Tri2Grid(X,Ver,step)

fii = -90:step:90;
lambdai = -180:step:180;

[fii,lambdai] = meshgrid(fii/180*pi, lambdai/180*pi);

vert1 = [X(1,Ver(1,:))' X(2,Ver(1,:))' X(3,Ver(1,:))'];
vert2 = [X(1,Ver(2,:))' X(2,Ver(2,:))' X(3,Ver(2,:))'];
vert3 = [X(1,Ver(3,:))' X(2,Ver(3,:))' X(3,Ver(3,:))'];


[xhat, yhat, zhat] = sph2cart(lambdai(:), fii(:),1);  
dir = [xhat, yhat, zhat]; 

orig = [0 0 0];
orig  = repmat(orig,size(vert1,1),1);

dirm = repmat(dir,fix(size(vert1,1)/size(dir,1)),1);

HowMuchMore = size(vert1)-size(dirm);

dirm = [dirm; dir(1:HowMuchMore,:)];

[intersect, t, u, v, xcoor] = TriangleRayIntersection (...
    orig, 10000000*dirm, vert1, vert2, vert3);


for i=1:numel(fii)
    
    [xhat, yhat, zhat] = sph2cart(lambdai(i), fii(i),1);   
    dir = [xhat, yhat, zhat];
    [intersect, t, u, v, xcoor] = TriangleRayIntersection (...
        orig, dir, vert1, vert2, vert3);
    progressbar(i/numel(fii));
      
end

progressbar(1);

x = xcoor(1,:);
y = xcoor(2,:);
z = xcoor(3,:);

[lambdai, fii, ri] = cart2sph(x,y,z);


