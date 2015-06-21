function V = Mesh2Volume(x,y,z)

s = size(x);
V = 0;

for i = 1:s(1)-1
    for j = 1:s(2)-1
        
        x1 = x(i,j);
        y1 = y(i,j);
        z1 = z(i,j);
        
        x2 = x(i+1,j);
        y2 = y(i+1,j);
        z2 = z(i+1,j);
        
        x3 = x(i,j+1);
        y3 = y(i,j+1);
        z3 = z(i,j+1);
        
        x4 = x(i+1,j+1);
        y4 = y(i+1,j+1);
        z4 = z(i+1,j+1);
        
        V1=(-x4.*y2.*z1+...
            x2.*y4.*z1+...
            x4.*y1.*z2-...
            x1.*y4.*z2-...
            x2.*y1.*z4+...
            x1.*y2.*z4)/6;
        
        V2=(-x3.*y2.*z1+...
            x2.*y3.*z1+...
            x3.*y1.*z2-...
            x1.*y3.*z2-...
            x2.*y1.*z3+...
            x1.*y2.*z3)/6;
        
        V = V + abs(V1) + abs(V2);
       
    end
end
