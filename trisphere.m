function [FV,P]=trisphere(level)
% TRISPHERE Generates a regular, triangular mesh of
%           the unit sphere surface. Starts with a
%           icosahedron, which is then tesselated.
%
%  CALL: [FV,P]=trisphere(level)
%
%        FV    = Face-vertex array, #{triangles} x 3
%        P     = Surface coordinates, as a 3 x #{vertices} matrix [X;Y;Z]
%        level = The number of tesselation iterations:
%                level 0 -> 20 triangles, 1->80, 2->320, 3->1280, etc.
%
% Example: [FV,P]=trisphere(3);trisurf(FV,P(1,:),P(2,:),P(3,:))

% Copyright Finn Lindgren 1997

% Build icosahedron.
P=zeros(3,12);
% Build top vertices.
P(:,1)=[0;0;1];
P(1,2:6)=cos(2*pi/5*(0:4))*2/sqrt(5);
P(2,2:6)=sin(2*pi/5*(0:4))*2/sqrt(5);
P(3,2:6)=ones(1,5)/sqrt(5);
% Build bottom vertices.
P(:,7:12)=[cos(pi/5) -sin(pi/5) 0;sin(pi/5) cos(pi/5) 0;0 0 -1]*P(:,1:6);

FV=[    ones(5,1)   [2:6]'    [3:6 2]'];
FV=[FV; ones(5,1)*7 12-[0:4]' 12-[1:4 0]'];
FV=[FV; 3 8 9; 4 9 10; 5 10 11; 6 11 12; 2 12 8];
FV=[FV; 8 3 2; 9 4 3; 10 5 4; 11 6 5; 12 2 6];

% Tesselate.
for lev=1:level
  [FV,P]=tessel(FV,P);
  l=sqrt(sum(P.^2));
  P=P./l(ones(3,1),:);
end


function [FV,P]=tessel(FV,P)
% TESSEL Divides all mesh triangles into 4 congruent subtriangles.
%
%  CALL: [FV,P]=tessel(FV,P)
%
%        FV,P : The mesh. (See TRISPHERE)

vv=vv2vv(fv2vv(FV));
vv=sort(vv')';
nP=size(P,2);
P=[P zeros(3,size(vv,1))];
for v=1:size(vv,1)
  P(:,nP+v)=(P(:,vv(v,1))+P(:,vv(v,2)))/2;
end

fv=FV;
nT=size(fv,1);
FV=zeros(nT*4,3);
for t=1:nT
  e=sort(fv(t,[1 2]));
  p1b=nP+find(vv(:,1)==e(1) & vv(:,2)==e(2));
  e=sort(fv(t,[2 3]));
  p2b=nP+find(vv(:,1)==e(1) & vv(:,2)==e(2));
  e=sort(fv(t,[3 1]));
  p3b=nP+find(vv(:,1)==e(1) & vv(:,2)==e(2));
  FV((t-1)*4+1,:)=[fv(t,1) p1b p3b];
  FV((t-1)*4+2,:)=[fv(t,2) p2b p1b];
  FV((t-1)*4+3,:)=[fv(t,3) p3b p2b];
  FV((t-1)*4+4,:)=[p1b p2b p3b];
end


function VV=fv2vv(FV,thetype,n)
% FV2VV Generates VV graph from FV graph.
%
%  CALL: VV = fv2vv(FV,1,n)
%        E  = fv2vv(FV,2), where E is a list of edges.

% Last modified 1997-06-06

if nargin<2
  thetype=1;
end
if thetype>2
  error('thetype must be 1 or 2. Did you enter  n  instead?')
end

if issparse(FV)
  if nargin<3
    n=size(FV,2);
  end
  FV=fv2fv(FV);
elseif nargin<3
  n=max(max(FV));
end

if isempty(FV)
  if isempty(n)
    n=0;
  end
  VV=sparse(n,n);
else
  VV=   sparse(FV(:,1),FV(:,2),ones(size(FV,1),1),n,n);
  VV=VV+sparse(FV(:,1),FV(:,3),ones(size(FV,1),1),n,n);
  VV=VV+sparse(FV(:,2),FV(:,3),ones(size(FV,1),1),n,n);
  VV=((VV+VV')>0);
end
if thetype==2
  VV=vv2vv(VV);
end


function FV=fv2fv(o_FV,n)
% FV2FV Translates between full/sparse representations of FV graph.
%
% FV=fv2fv(o_FV),   o_FV sparse
% FV=fv2fv(o_FV,n), o_FV full

% Last modified 1997-06-06

if issparse(o_FV)
  FV=zeros(size(o_FV,1),sum(o_FV(1,:)>0));
  for t=1:size(o_FV,1)
    i=find(o_FV(t,:))';
    FV(t,o_FV(t,i))=i;
  end
else
  if isempty(o_FV)
    FV=[];
  else
    if nargin<2
      n=max(max(o_FV));
    end
    FV=sparse(size(o_FV,1),n);
    for t=1:size(o_FV,1)
      FV(t,o_FV(t,:))=1:size(o_FV,2);
    end
  end
end


function VV=vv2vv(o_VV,n)
% VV2VV Translates between full/sparse representations of VV graph.
%
% VV = vv2vv(E,n)
% E  = vv2vv(VV)

if issparse(o_VV)
  [i,j]=find(triu(o_VV));
  VV=[i j];
  return;
end

if nargin<2
  n=max(o_VV(:));
end
VV=sparse(o_VV(:,1),o_VV(:,2),ones(size(o_VV,1),1),n,n);
VV=(VV+VV')>0;