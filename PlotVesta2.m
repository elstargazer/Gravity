filename='Vesta20130522shape_geoc_elev_12m_gridline.grd';
[lambda,fi,r]=ReadGRD(filename);

[x,y,z]=sph2cart(lambda/180*pi,fi/180*pi,r);


at=280.86977;
bt=226.49124;
et=Eccentricity(at,bt);

[~,~,H]=XYZ2BLH(x,y,z,at,et);

AGUaxes;

surfm(fi/180*pi,lambda/180*pi,H);
