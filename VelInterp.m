function v=VelInterp(a1,a2,v1,v2,a1i,a2i)

v1i=griddata(a1,a2,v1,a1i,a2i,'cubic');
v2i=griddata(a1,a2,v2,a1i,a2i,'cubic');

v=sqrt(v1i.^2+v2i.^2);



