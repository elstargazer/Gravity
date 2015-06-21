% function SplitInThree(FV,ver)

zmax = max(ver(:,3));
zmin = min(ver(:,3));

range = zmax - zmin;

zS = zmin+range/9;
zN = zmax-range/16;

sFV = size(FV);

s1=0;
s2=0;
s3=0;

FV1=zeros(size(FV));
FV2=FV1;
FV3=FV1;

i1=1;
i2=1;
i3=1;

for i=1:sFV(1)
    
    cond1 = (ver(FV(i,1),3) > zN) && ...
            (ver(FV(i,2),3) > zN) && ...
            (ver(FV(i,3),3) > zN);
        
    cond2 = (ver(FV(i,1),3) < zS) && ...
            (ver(FV(i,2),3) < zS) && ...
            (ver(FV(i,3),3) < zS);
        
    cond3 = (~cond1) && (~cond2);
    
    s1=s1+cond1;
    s2=s2+cond2;
    s3=s3+cond3;
    
    if cond1
        FV1(i1,:) = FV(i,:);
        i1=i1+1;
    end
    
    if cond2
        FV2(i2,:) = FV(i,:);
        i2=i2+1;
    end
    
    if cond3
        FV3(i3,:) = FV(i,:);
        i3=i3+1;
    end
    
    
end

FV1=FV1(1:s1,:);
FV2=FV2(1:s2,:);
FV3=FV3(1:s3,:);


figure; hold on;
trisurf(FV1,ver(:,1),ver(:,2),ver(:,3),r); shading interp
trisurf(FV2,ver(:,1),ver(:,2),ver(:,3),r); shading interp
trisurf(FV3,ver(:,1),ver(:,2),ver(:,3),r); shading interp
colorbar; axis equal

A1.faces=FV1;
A1.vertices=ver;

A2.faces=FV2;
A2.vertices=ver;

A3.faces=FV3;
A3.vertices=ver;

filename_bin1=['VestaShapeModel_Dinara_final2_1_' num2str(NTess) '_' ...
     '3m' '_bin.stl'];
filename_bin2=['VestaShapeModel_Dinara_final2_2_' num2str(NTess) '_' ...
     '3m' '_bin.stl'];
filename_bin3=['VestaShapeModel_Dinara_final2_3_' num2str(NTess) '_' ...
     '3m' '_bin.stl'];

stlwrite(filename_bin1, A1,'MODE','binary');
stlwrite(filename_bin2, A2,'MODE','binary');
stlwrite(filename_bin3, A3,'MODE','binary');







