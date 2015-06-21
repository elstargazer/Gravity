ShapeModelFile='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_SHAPE_POST_VESTA/SHAPE.TXT';

data=load(ShapeModelFile);

x=data(:,1);
y=data(:,2);
z=data(:,3);

P=[x'; y'; z'];

K = convhulln(P');  
K = unique(K(:));  
Q = P(:,K);
[A, c] = MinVolEllipse(Q, .01)


[U Q V] = svd(A);

a = 1/sqrt(Q(1,1)); 
b = 1/sqrt(Q(2,2)); 
c = 1/sqrt(Q(3,3));

Volume_min_ell=4/3*pi*a*b*c;

R_b=292.800;

Volume_min_sph=4/3*pi*(R_b)^3;

Volume_min_ell/Volume_min_sph

