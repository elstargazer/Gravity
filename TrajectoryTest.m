% ccc


%% Read Sumfiles

filelist='SumfilesList';
folder='/Users/antonermakov/Dawn/Vesta_Shape_Models/Gaskell/GASKELL_CLAUDIA_2013_05_22/GASKELL_CLAUDIA_2013_05_22_SUMFILES';

[x,y,z,t]=ReadSumfileList(folder,filelist);

r=sqrt(x.*x+y.*y+z.*z);


%% Load spice kernels

metakernel='/Users/antonermakov/Dawn/SPICE/meta/metakernel_rec.txt';

cspice_furnsh(metakernel);


%% Gettign coordinates

mu=17.28;

target   = '-203';
epoch0    = 'June 1, 2011 12:00 AM PST';
epoch1    = 'March 1, 2012 12:00 AM PST';
frame    = 'J2000';
abcorr   = 'none';
observer = '2000004';

et0 = cspice_str2et( epoch0 );
et1 = cspice_str2et( epoch1 );

et_step=60*5;

et=et0:et_step:et1;
                                 
[ state, ltime ] = cspice_spkezr( target, et, frame, ...
                                   abcorr, observer);      
                                                            
[a,e,i,w,W,M]=koor_el2(state(1,:),state(2,:),state(3,:),...
                      state(4,:),state(5,:),state(6,:),mu);                          
                                 
% plot orbit                             
figure; hold on;

plot3(state(1,:),state(2,:),state(3,:),'b-');

[xv,yv,zv]=sphere(100);
Rv=261;
xv=xv*Rv;
yv=yv*Rv;
zv=zv*Rv;

surf(xv,yv,zv);


%% Compare

target   = '2000004';
epoch0    = 'June 1, 2011 12:00 AM PST';
epoch1    = 'March 1, 2012 12:00 AM PST';
frame    = 'IAU_VESTA';
abcorr   = 'none';
observer = '-203';

et_image=cspice_str2et( t);

[ pos_image, ltime ] = cspice_spkpos( target, et_image, frame, ...
                                   abcorr, observer);   

% plot orbit                               
                    
figure; hold on;

plot3(state(1,:),state(2,:),state(3,:),'b.','MarkerSize',10);


dx=x'-pos_image(1,:);
dy=y'-pos_image(2,:);
dz=z'-pos_image(3,:);

dr=sqrt(dx.*dx+dy.*dy+dz.*dz);


figure; hold on;
plot(dx,'r.','MarkerSize',10);

figure; hold on;
plot(dy,'r.','MarkerSize',10);

figure; hold on;
plot(dz,'r.','MarkerSize',10);

figure; hold on;
plot(dr,'b.','MarkerSize',10);

set(gca,'YScale','log');













                 
                               
