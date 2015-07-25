ccc

filename = 'mesh.txt';


x = -1:0.25:1;
y = -1:0.25:1;
z = -1:0.25:1;

nelems_row = numel(x)-1;

d = 1;

x=x*d;
y=y*d;
z=z*d;

nelems = (numel(x)-1).^3;
con = zeros(nelems,8);
[x,y,z] = ndgrid(x,y,z);

for i = 1:nelems
    con(i,:) = [1 2 3 4 5 6 7 8];
end

figure; hold on; box on;
plot3(x(:),y(:),z(:),'o','MarkerSize',10,'MarkerFaceColor','k');
axis equal;
xlim(1.1*[-1 1]);
ylim(1.1*[-1 1]);
zlim(1.1*[-1 1]);


i=2;
j=2;
k=2;

[ x(i,j,k) y(i,j,k) z(i,j,k)]+1

% 1 - 1 1 1
% 4 - 1 1 2
% 5 - 1 2 1
% 8 - 1 2 2
% 2 - 2 1 1
% 3 - 2 1 2
% 6 - 2 2 1
% 7 - 2 2 2

% 1-2
i1s = 0; i1f = 1;
j1s = 0; j1f = 0;
k1s = 0; k1f = 0;

% 2-3
i2s = 1; i2f = 1;
j2s = 0; j2f = 0;
k2s = 0; k2f = 1;

% 3-4
i3s = 1; i3f = 0;
j3s = 0; j3f = 0;
k3s = 1; k3f = 1;

% 4-1
i4s = 0; i4f = 0;
j4s = 0; j4f = 0;
k4s = 1; k4f = 0;

% 1-5
i5s = 0; i5f = 0;
j5s = 0; j5f = 1;
k5s = 0; k5f = 0;

% 2-6
i6s = 1; i6f = 1;
j6s = 0; j6f = 1;
k6s = 0; k6f = 0;

% 3-7
i7s = 1; i7f = 1;
j7s = 0; j7f = 1;
k7s = 1; k7f = 1;

% 4-8
i8s = 0; i8f = 0;
j8s = 0; j8f = 1;
k8s = 1; k8f = 1;

% 5-6
i9s = 0; i9f = 1;
j9s = 1; j9f = 1;
k9s = 0; k9f = 0;

% 6-7
i10s = 1; i10f = 1;
j10s = 1; j10f = 1;
k10s = 0; k10f = 1;

% 7-8
i11s = 1; i11f = 0;
j11s = 1; j11f = 1;
k11s = 1; k11f = 1;

% 8-5
i12s = 0; i12f = 0;
j12s = 1; j12f = 1;
k12s = 1; k12f = 0;

for i=1:nelems_row
    for j=1:nelems_row
        for k=1:nelems_row
            
            line([x(i+i1s,j+j1s,k+k1s) x(i+i1f,j+j1f,k+k1f)],...
                [y(i+i1s,j+j1s,k+k1s) y(i+i1f,j+j1f,k+k1f)],...
                [z(i+i1s,j+j1s,k+k1s) z(i+i1f,j+j1f,k+k1f)]);
            
            line([x(i+i2s,j+j2s,k+k2s) x(i+i2f,j+j2f,k+k2f)],...
                [y(i+i2s,j+j2s,k+k2s) y(i+i2f,j+j2f,k+k2f)],...
                [z(i+i2s,j+j2s,k+k2s) z(i+i2f,j+j2f,k+k2f)]);
            
            line([x(i+i3s,j+j3s,k+k3s) x(i+i3f,j+j3f,k+k3f)],...
                [y(i+i3s,j+j3s,k+k3s) y(i+i3f,j+j3f,k+k3f)],...
                [z(i+i3s,j+j3s,k+k3s) z(i+i3f,j+j3f,k+k3f)]);
            
            line([x(i+i4s,j+j4s,k+k4s) x(i+i4f,j+j4f,k+k4f)],...
                [y(i+i4s,j+j4s,k+k4s) y(i+i4f,j+j4f,k+k4f)],...
                [z(i+i4s,j+j4s,k+k4s) z(i+i4f,j+j4f,k+k4f)]);
            
            line([x(i+i5s,j+j5s,k+k5s) x(i+i5f,j+j5f,k+k5f)],...
                [y(i+i5s,j+j5s,k+k5s) y(i+i5f,j+j5f,k+k5f)],...
                [z(i+i5s,j+j5s,k+k5s) z(i+i5f,j+j5f,k+k5f)]);
            
            line([x(i+i6s,j+j6s,k+k6s) x(i+i6f,j+j6f,k+k6f)],...
                [y(i+i6s,j+j6s,k+k6s) y(i+i6f,j+j6f,k+k6f)],...
                [z(i+i6s,j+j6s,k+k6s) z(i+i6f,j+j6f,k+k6f)]);
            
            line([x(i+i7s,j+j7s,k+k7s) x(i+i7f,j+j7f,k+k7f)],...
                [y(i+i7s,j+j7s,k+k7s) y(i+i7f,j+j7f,k+k7f)],...
                [z(i+i7s,j+j7s,k+k7s) z(i+i7f,j+j7f,k+k7f)]);
            
            line([x(i+i8s,j+j8s,k+k8s) x(i+i8f,j+j8f,k+k8f)],...
                [y(i+i8s,j+j8s,k+k8s) y(i+i8f,j+j8f,k+k8f)],...
                [z(i+i8s,j+j8s,k+k8s) z(i+i8f,j+j8f,k+k8f)]);
            
            line([x(i+i9s,j+j9s,k+k9s) x(i+i9f,j+j9f,k+k9f)],...
                [y(i+i9s,j+j9s,k+k9s) y(i+i9f,j+j9f,k+k9f)],...
                [z(i+i9s,j+j9s,k+k9s) z(i+i9f,j+j9f,k+k9f)]);
            
            line([x(i+i10s,j+j10s,k+k10s) x(i+i10f,j+j10f,k+k10f)],...
                [y(i+i10s,j+j10s,k+k10s) y(i+i10f,j+j10f,k+k10f)],...
                [z(i+i10s,j+j10s,k+k10s) z(i+i10f,j+j10f,k+k10f)]);
            
            line([x(i+i11s,j+j11s,k+k11s) x(i+i11f,j+j11f,k+k11f)],...
                [y(i+i11s,j+j11s,k+k11s) y(i+i11f,j+j11f,k+k11f)],...
                [z(i+i11s,j+j11s,k+k11s) z(i+i11f,j+j11f,k+k11f)]);
            
            line([x(i+i12s,j+j12s,k+k12s) x(i+i12f,j+j12f,k+k12f)],...
                [y(i+i12s,j+j12s,k+k12s) y(i+i12f,j+j12f,k+k12f)],...
                [z(i+i12s,j+j12s,k+k12s) z(i+i12f,j+j12f,k+k12f)]);
        end
    end
end



%%

xs = x .* sqrt(1.0 - (y.*y/2.0) - (z.*z/2.0) + (y.*y.*z.*z/3.0));
ys = y .* sqrt(1.0 - (z.*z/2.0) - (x.*x/2.0) + (z.*z.*x.*x/3.0));
zs = z .* sqrt(1.0 - (x.*x/2.0) - (y.*y/2.0) + (x.*x.*y.*y/3.0));

rs = sqrt(xs.*xs + ys.*ys + zs.*zs);

x = xs;
y = ys;
z = zs;

% figure;
% plot(rs(:));

figure; hold on;box on
plot3(xs(:),ys(:),zs(:),'o','MarkerSize',10,'MarkerFaceColor','k');
axis equal;

xlim(1.1*[-1 1]);
ylim(1.1*[-1 1]);
zlim(1.1*[-1 1]);



for i=1:nelems_row
    for j=1:nelems_row
        for k=1:nelems_row
            
            line([x(i+i1s,j+j1s,k+k1s) x(i+i1f,j+j1f,k+k1f)],...
                [y(i+i1s,j+j1s,k+k1s) y(i+i1f,j+j1f,k+k1f)],...
                [z(i+i1s,j+j1s,k+k1s) z(i+i1f,j+j1f,k+k1f)]);
            
            line([x(i+i2s,j+j2s,k+k2s) x(i+i2f,j+j2f,k+k2f)],...
                [y(i+i2s,j+j2s,k+k2s) y(i+i2f,j+j2f,k+k2f)],...
                [z(i+i2s,j+j2s,k+k2s) z(i+i2f,j+j2f,k+k2f)]);
            
            line([x(i+i3s,j+j3s,k+k3s) x(i+i3f,j+j3f,k+k3f)],...
                [y(i+i3s,j+j3s,k+k3s) y(i+i3f,j+j3f,k+k3f)],...
                [z(i+i3s,j+j3s,k+k3s) z(i+i3f,j+j3f,k+k3f)]);
            
            line([x(i+i4s,j+j4s,k+k4s) x(i+i4f,j+j4f,k+k4f)],...
                [y(i+i4s,j+j4s,k+k4s) y(i+i4f,j+j4f,k+k4f)],...
                [z(i+i4s,j+j4s,k+k4s) z(i+i4f,j+j4f,k+k4f)]);
            
            line([x(i+i5s,j+j5s,k+k5s) x(i+i5f,j+j5f,k+k5f)],...
                [y(i+i5s,j+j5s,k+k5s) y(i+i5f,j+j5f,k+k5f)],...
                [z(i+i5s,j+j5s,k+k5s) z(i+i5f,j+j5f,k+k5f)]);
            
            line([x(i+i6s,j+j6s,k+k6s) x(i+i6f,j+j6f,k+k6f)],...
                [y(i+i6s,j+j6s,k+k6s) y(i+i6f,j+j6f,k+k6f)],...
                [z(i+i6s,j+j6s,k+k6s) z(i+i6f,j+j6f,k+k6f)]);
            
            line([x(i+i7s,j+j7s,k+k7s) x(i+i7f,j+j7f,k+k7f)],...
                [y(i+i7s,j+j7s,k+k7s) y(i+i7f,j+j7f,k+k7f)],...
                [z(i+i7s,j+j7s,k+k7s) z(i+i7f,j+j7f,k+k7f)]);
            
            line([x(i+i8s,j+j8s,k+k8s) x(i+i8f,j+j8f,k+k8f)],...
                [y(i+i8s,j+j8s,k+k8s) y(i+i8f,j+j8f,k+k8f)],...
                [z(i+i8s,j+j8s,k+k8s) z(i+i8f,j+j8f,k+k8f)]);
            
            line([x(i+i9s,j+j9s,k+k9s) x(i+i9f,j+j9f,k+k9f)],...
                [y(i+i9s,j+j9s,k+k9s) y(i+i9f,j+j9f,k+k9f)],...
                [z(i+i9s,j+j9s,k+k9s) z(i+i9f,j+j9f,k+k9f)]);
            
            line([x(i+i10s,j+j10s,k+k10s) x(i+i10f,j+j10f,k+k10f)],...
                [y(i+i10s,j+j10s,k+k10s) y(i+i10f,j+j10f,k+k10f)],...
                [z(i+i10s,j+j10s,k+k10s) z(i+i10f,j+j10f,k+k10f)]);
            
            line([x(i+i11s,j+j11s,k+k11s) x(i+i11f,j+j11f,k+k11f)],...
                [y(i+i11s,j+j11s,k+k11s) y(i+i11f,j+j11f,k+k11f)],...
                [z(i+i11s,j+j11s,k+k11s) z(i+i11f,j+j11f,k+k11f)]);
            
            line([x(i+i12s,j+j12s,k+k12s) x(i+i12f,j+j12f,k+k12f)],...
                [y(i+i12s,j+j12s,k+k12s) y(i+i12f,j+j12f,k+k12f)],...
                [z(i+i12s,j+j12s,k+k12s) z(i+i12f,j+j12f,k+k12f)]);
        end
    end
end

%% Write file

in = fopen(filename,'w');

fprintf(in,'%23.16E, %23.16E, %23.16E\n',x(:),y(:),z(:));

s = size(x);
ver=zeros(1,8);

for i=1:nelems_row
    for j=1:nelems_row
        for k=1:nelems_row
    
            ver(1) = sub2ind(s, i, j, k);
            ver(2) = sub2ind(s, i+1, j, k);
            ver(3) = sub2ind(s, i+1, j, k+1);
            ver(4) = sub2ind(s, i, j, k+1);
            ver(5) = sub2ind(s, i, j+1, k);
            ver(6) = sub2ind(s, i+1, j+1, k);
            ver(7) = sub2ind(s, i+1, j+1, k+1);
            ver(8) = sub2ind(s, i, j+1, k+1);
    
            fprintf(in,'%d, %d, %d, %d, %d, %d, %d, %d\n',ver);        
        end
    end
end

fclose(in);





