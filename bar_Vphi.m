clear;
clear;
ffname='RealPlutinosNpl';
fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

fontsize=15;
N=50;
yylim=0.1;

set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
%     r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/r2hill_record.txt'));
%     r2hill_mean=mean(r2hill_record);
%     hill=sqrt(r2hill_mean);
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_tp.txt'));
%     AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_pl.txt'));
%     AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_tp.txt'));

    %di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/di_record_inout.txt'));
    Vp=PV_record_pl(:,4:6);
    Vt=PV_record_tp(:,4:6);
    Vpnorm=(Vp(:,1).^2+Vp(:,2).^2+Vp(:,3).^2).^(1/2);
    Vtnorm=(Vt(:,1).^2+Vt(:,2).^2+Vt(:,3).^2).^(1/2);
    dot=Vp(:,1).*Vt(:,1)+Vp(:,2).*Vt(:,2)+Vp(:,3).*Vt(:,3);
    cosphi=dot./(Vpnorm.*Vtnorm);
    cosphix=-1:2/N:1;
    cosphix=cosphix';
    countx=histcounts(cosphi,cosphix)'/length(Vp);
    subplot(2,2,isub);
    plot(cosphix(2:end),countx,'k.-');
end