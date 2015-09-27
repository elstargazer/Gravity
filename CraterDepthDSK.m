function [depth,depth_std,rc_e]=CraterDepthDSK(shape_filename,latcc,loncc,rc)


%% Setup shape
TOL      =  1.d-12;

% Open the DSK file for read access.
% We use the DAS-level interface for
% this function.
%
handle = cspice_dasopr( shape_filename );

%
% Begin a forward search through the
% kernel, treating the file as a DLA.
% In this example, it's a very short
% search.
%
[dladsc, found] = cspice_dlabfs( handle );

if ~found
    
    %
    % We arrive here only if the kernel
    % contains no segments. This is
    % unexpected, but we're prepared for it.
    %
    fprintf( 'No segments found in DSK file %s\n', dsk )
    return
    
end

%
% If we made it this far, DLADSC is the
% DLA descriptor of the first segment.


%% Make grid
step = 2;
L = 20;
lat=(-90:step:90) * cspice_rpd();
lon=(0:step:360) * cspice_rpd();

[lati, loni] = meshgrid(lat, lon);

s = size( lati );

lati = lati(:);
loni = loni(:);

grid = zeros(2,numel(lati));

grid(1,:) = loni;
grid(2,:) = lati;

%
[spoints, ~] = cspice_llgrid_pl02( handle, dladsc, grid );

x = spoints(1,:);
y = spoints(2,:);
z = spoints(3,:);

x = reshape(x,s);
y = reshape(y,s);
z = reshape(z,s);

lati = reshape(lati,s);
loni = reshape(loni,s);

ri = sqrt(x.*x+y.*y+z.*z);

lmcosi = xyz2plm(ri',L);

ri_sh = plm2xyz(lmcosi,step);
ri_sh = ri_sh';


%% Setup parameters

cf1=1.4;
cf2=1.3; % cf2 should be less than cf1

R=lmcosi(1,3);
Nprof=30;
npts=100;

botloc=0.25;
aroundrim=0.1;

Ncr=numel(latcc);
depth=NaN(1,Ncr);
depth_std=NaN(1,Ncr);
rc_e=NaN(1,Ncr);

progressbar(0);

az=linspace(0,180,Nprof);

for i=1:Ncr    
    clear lattrk lontrk 
    for j=1:Nprof
        
        [latpf,lonpf] = reckon(latcc(i),loncc(i),rc(i)*cf2,az(j));
        [latpb,lonpb] = reckon(latcc(i),loncc(i),rc(i)*cf2,180+az(j));
        [lattrk(:,j),lontrk(:,j)] = track([latpf latpb],[lonpf lonpb],[],...
            'degrees',npts);
        %         plotm(lattrk/180*pi,lontrk/180*pi,'-','LineWidth',1)
    end
    
    lattrk(end,:) = [];
    lontrk(end,:) = [];
     
    lonrim=[];
    latrim=[];
    
    lonbot=[];
    latbot=[];
    
    rim_peaks=[];
    crat_bot=[];
    
    grid_prof = zeros(2,numel(lattrk));
    
    s_prof = size(lattrk);
    grid_prof(1,:) = lontrk(:);
    grid_prof(2,:) = lattrk(:);
    
    [spoints_prof, ~] = cspice_llgrid_pl02( handle, dladsc, grid_prof );
    
    x_prof = spoints_prof(1,:);
    y_prof = spoints_prof(2,:);
    z_prof = spoints_prof(3,:);
    
    r_prof = reshape(sqrt(x_prof.*x_prof+y_prof.*y_prof+z_prof.*z_prof),s_prof);  
    r_prof_sh = reshape(plm2xyz(lmcosi,lattrk(:),lontrk(:)),s_prof);

    Ht = r_prof - r_prof_sh;
     
    for j=1:Nprof
        
        % finding rims
        [pks,locs] = findpeaks(Ht(:,j));
        %         plot(locs,pks,'x','MarkerSize',3);
        % exclude peaks in the center
        
        %         exccen=(abs(locs-(npts/2))<(npts*botloc/cf2));
   
        r1=fix(npts*(cf2-1-aroundrim)/(2*cf2));
        r2=fix(npts*(cf2-1+aroundrim)/(2*cf2))+1;
        
        r3=fix(npts*(cf2+1-aroundrim)/(2*cf2));
        r4=fix(npts*(cf2+1+aroundrim)/(2*cf2))+1;
        
        rimonly=((locs>r1) & (locs<r2)) | ((locs>r3) & (locs<r4));
        
        pks=pks(rimonly);
        locs=locs(rimonly);
        
        %         figure(f4)
        %         plot(locs,pks,'^','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','w');
               
        lonrim=[lonrim; lontrk(locs,j)];
        latrim=[latrim; lattrk(locs,j)];
        rim_peaks=[rim_peaks; pks];
        
        
        % finding bottom
        [bots,botin]=min(Ht(:,j));
        
        %         figure(f4)
        %         plot(botin,bots,'v','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','w');
        
        latbot=[latbot lattrk(botin,j)];
        lonbot=[lonbot lontrk(botin,j)];
        
        crat_bot=[crat_bot bots];
        
    end
    
    [xsrim,ysrim]=sph2stereo(lonrim/180*pi,latrim/180*pi,loncc(i)/180*pi,latcc(i)/180*pi,R);
    [xsbot,ysbot]=sph2stereo(lonbot/180*pi,latbot/180*pi,loncc(i)/180*pi,latcc(i)/180*pi,R);
    
%     figure(f3)
%     plot3(xsrim,ysrim,rim_peaks+350,'^','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','w');
%     plot3(xsbot,ysbot,crat_bot+350,'v','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','w');
%     
%     cbar=colorbar('FontSize',20,'Color','w');
%     ylabel(cbar,'Elevation [m]','FontSize',20,'Color','w');
%     
%     plot3m(latrim/180*pi, lonrim/180*pi,3,'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','r');
%     
    % depth
    try
        %         b=min(Ht(:));
        b=median(crat_bot);
        t=median(rim_peaks);
        
        b_std=std(crat_bot);
        t_std=std(rim_peaks);
        
        depth(i)=t-b;
        depth_std(i)=sqrt(t_std.*t_std+b_std.*b_std);
    catch
        depth(i)=NaN;
        depth_std(i)=NaN;
        
    end
    
    % radius equivalent
    try
        ellipse_t=fit_ellipse(xsrim,ysrim);
        ac=ellipse_t.long_axis;
        cc=ellipse_t.short_axis;
        % Htm=mean(Ht,2);
        rc_e(i)=sqrt(ac(i)*cc(i));
    catch
        rc_e(i)=NaN;
    end
    
    progressbar(i/Ncr);
    
end

progressbar(1);
