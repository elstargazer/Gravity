ccc

ell_thomas=[487.3 487.3 454.7];
ell_carry=[479.7 479.7 444.4];
ell_occ=[479.6 479.6 453.4];
ell_occ2=[481.6 481.6 450.1];
ell_dlr=[483.26 480.97 446.13];
ell_jpl=[484.17 481.41 447.80];
ell_limb=[484.5 482.5 447.6];

ell_thomas_err=[1.8 1.8 1.6];
ell_carry_err=[2.3 2.3 2.1];
ell_occ_err=[2.4 2.4 4.5];
ell_occ_err2=[2.4 2.4 2.0];
ell_dlr_err=0*[2.3 2.3 2.1]+0.1;
ell_jpl_err=0*[1.8 1.8 1.6]+0.1;
ell_limb_err=[0.7 0.7 0.4];

rvol_thomas=prod(ell_thomas).^(1/3);
rvol_carry=prod(ell_carry).^(1/3);
rvol_occ=prod(ell_occ).^(1/3);
rvol_occ2=prod(ell_occ2).^(1/3);
rvol_jpl=prod(ell_jpl).^(1/3);
rvol_dlr=prod(ell_dlr).^(1/3);
rvol_limb=prod(ell_limb).^(1/3);

[fp_thomas,~,~,fp_err_thomas]=abc2f(ell_thomas,ell_thomas_err);
[fp_carry,~,~,fp_err_carry]=abc2f(ell_carry,ell_carry_err);
[fp_occ,~,~,fp_err_occ]=abc2f(ell_occ,ell_occ_err);
[fp_occ2,~,~,fp_err_occ2]=abc2f(ell_occ2,ell_occ_err2);



rvol_thomas_err=1.3;
rvol_carry_err=1.6;
rvol_occ_err=1.8;
rvol_occ_err2=1.3;
rvol_jpl_err=.1;
rvol_dlr_err=.1;
rvol_limb_err=.3;

l=[1 2 3];

figure; hold on; 
set(gca,'FontSize',20);
box on;

plot(l-0.2,ell_thomas,'or','MarkerSize',5);
plot(l-0.1,ell_carry,'ob','MarkerSize',5);
plot(l+0.0,ell_occ,'oy','MarkerSize',5);
plot(l+0.1,ell_occ2,'og','MarkerSize',5);
% plot(l+0.2,ell_dlr,'oc','MarkerSize',5);
% plot(l+0.3,ell_jpl,'om','MarkerSize',5);
% plot(l+0.4,ell_limb,'ok','MarkerSize',5);


errorbar(l-0.2,ell_thomas,ell_thomas_err,'.r','LineWidth',3);
errorbar(l-0.1,ell_carry,ell_carry_err,'.b','LineWidth',3);
errorbar(l+0.0,ell_occ,ell_occ_err,'.y','LineWidth',3);
errorbar(l+0.1,ell_occ2,ell_occ_err2,'.g','LineWidth',3);
% errorbar(l+0.2,ell_dlr,ell_dlr_err,'.c','LineWidth',3);
% errorbar(l+0.3,ell_jpl,ell_jpl_err,'.m','LineWidth',3);
% errorbar(l+0.4,ell_limb,ell_limb_err,'.k','LineWidth',3);

%
plot(4-0.2,rvol_thomas,'or','MarkerSize',5);
plot(4-0.1,rvol_carry,'ob','MarkerSize',5);
plot(4+0.0,rvol_occ,'oy','MarkerSize',5);
plot(4+0.1,rvol_occ2,'og','MarkerSize',5);
% plot(4+0.2,rvol_dlr,'oc','MarkerSize',5);
% plot(4+0.3,rvol_jpl,'om','MarkerSize',5);
% plot(4+0.4,rvol_limb,'ok','MarkerSize',5);


errorbar(4-0.2,rvol_thomas,rvol_thomas_err,'.r','LineWidth',3);
errorbar(4-0.1,rvol_carry,rvol_carry_err,'.b','LineWidth',3);
errorbar(4+0.0,rvol_occ,rvol_occ_err,'.y','LineWidth',3);
errorbar(4+0.1,rvol_occ2,rvol_occ_err2,'.g','LineWidth',3);
% errorbar(4+0.2,rvol_dlr,rvol_dlr_err,'.c','LineWidth',3);
% errorbar(4+0.3,rvol_jpl,rvol_jpl_err,'.m','LineWidth',3);
% errorbar(4+0.4,rvol_limb,rvol_limb_err,'.k','LineWidth',3);

% legend({'Thomas et al., 2005 (HST)','Carry et al., 2008 (Keck)','Millis et al., 1987 (Occ)','Millis et al., 1987 (Occ2)','Dawn SPG (DLR)','Dawn SPC (JPL)','Dawn limb (DLR)'},'FontSize',20);

legend({'Thomas et al., 2005 (HST)',...
    'Carry et al., 2008 (Keck)',...
    'Millis et al., 1987 (Occ, rmin)',...
    'Millis et al., 1987 (Occ, tmin)'...
    },'FontSize',20);

ylabel('Length [km]','FontSize',20);

grid on;

xlim([0 5]);

set(gca,'XTickLabel',{' ','a','b','c','r_{eqvol}',' '});