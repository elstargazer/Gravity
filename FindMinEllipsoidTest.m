ccc

filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
NTess=5;

[TR,x1,y1,z1]=Grid2TriShapeModel(filename,NTess);

% Brilloiun sphere

rb=293.000;
Vshb=4/3*pi*rb^3;

% Park et al 2013 Ellipsoid
ap=303.860;
bp=288.660;
cp=246.460;

Vp=4/3*pi*ap*bp*cp;

Tol=1e-2;

P=[x1'; y1'; z1'];
K = convhulln(P');  
K = unique(K(:));  
Qc = P(:,K);

[A , c] = MinVolEllipse(Qc, Tol);

[U, Q, V] = svd(A);

semiaxes=1./sqrt(diag(Q));

Va=4/3*pi*prod(semiaxes)

