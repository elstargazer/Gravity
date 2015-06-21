function WriteInputTopoSHFile(lmcosi,MaxDegree,GM)
% MaxDegree=50;
Rref=lmcosi(1,3);

lmcosi(:,3)=lmcosi(:,3)/Rref;
lmcosi(:,4)=lmcosi(:,4)/Rref;

invf=1;
% GM=17.28e9;
omega=1;

NumberOfCoef=numel(lmcosi(:,1));

InputFileName='/Users/antonermakov/Dawn/Balmino/VestaTest/file_in_harmo_shape_topograv';

in=fopen(InputFileName,'w');

Comment=['VESTA TOPO HARMONIC COEF. (MAXDEGREE=' num2str(MaxDegree) ')'];
Title1='         ae                  1/f                 gm                 omega';
Papameters=[num2str(Rref,'%1.14E') num2str(invf,'%1.14E') num2str(GM,'%1.14E') num2str(omega,'%1.14E')];
Info='reference date : 2012.00';
Title2=['maximal degree :' num2str(MaxDegree,'%62i') '     sigmas calibration factor : .1000e+01 (already applied)'];
Title3=' l  m dot         cbar                sbar             sigma c      sigma s  lib';
lib=0;

fprintf(in,'%s\n',Comment);
fprintf(in,'%s\n',Title1);
fprintf(in,'%s\n',Papameters);
fprintf(in,'%s\n',Info);
fprintf(in,'%s\n',Title2);
fprintf(in,'%s\n',Title3);

sigmac=lmcosi(:,3)*0;
sigmas=lmcosi(:,4)*0;
lib=lmcosi(:,4)*0;

Degree=[];
Order=[];

for i=0:MaxDegree    
   Degree=[Degree i:MaxDegree];    
   Order=[Order i*ones(1,MaxDegree-i+1)];    
end
Degree=Degree';
Order=Order';


IndexDegree=find(lmcosi(:,1)==Degree);
IndexOrder=find(lmcosi(:,2)==Order);

for i=1:NumberOfCoef    

    Index=(  (lmcosi(:,1)==Degree(i)) & (lmcosi(:,2)==Order(i))  );    
    Degree(i);
    Order(i);
    find(Index)    
    lmcosi_r(i,:)=[Degree(i) Order(i) lmcosi(Index,3) lmcosi(Index,4)];    
    
end

lmcosi_r(1:MaxDegree+1,4)=0;

fprintf(in,'%3i%3i   %+1.14E%+1.14E%+1.6E%+1.6E  %1i\n',[lmcosi_r(:,1),lmcosi_r(:,2),lmcosi_r(:,3),lmcosi_r(:,4),sigmac,sigmas,lib]');




