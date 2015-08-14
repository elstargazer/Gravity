function [VT,ET,FTb,FT,faceBoundaryMarker,cmT] = MergeNodeSets(V2,V1,E2,E1,cm2,cm1,F1)

% Merge two meshes from inside out
% V1 and E1 are vertices and elements of the outer layer
% V2 and E2 are vertices and elements of the inner layer
% F1 is boundary of the outer layer


% Merging node sets
VT=[V2;V1];
ET=[E2;E1+size(V2,1)];
cmT = [cm2; cm1];
[~,ind1,ind2]=unique(pround(VT,5),'rows');
VT=VT(ind1,:);
ET=ind2(ET);

FTb=ind2(F1+size(V2,1));
faceBoundaryMarker=ones(size(F1,1),1);

%Get faces
[FT,~]=element2patch(ET,[],'hex8');