%% CE location distribution
clear;
ffname1='RealPlutinosNpl';
ffname2='RealTrojans';
fname11='1999CE119_1Gyr_40pl';
fname12='2001FU172_1Gyr';
%fname='2000CK105_1Gyr';
%fname='2004KV18_1Gyr';

fname21='2006RJ103_1Gyr';
fname22='2001FU172&2006RJ103_1Gyr';


% 
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname12,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname12,'/PV_record_tp.txt'));
% PV_rlt12=PV_record_tp-PV_record_pl;
% r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname12,'/r2hill_record.txt'));
% r2hill_mean=mean(r2hill_record);
% hill12=sqrt(r2hill_mean);
% 
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname21,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname21,'/PV_record_tp.txt'));
% PV_rlt21=PV_record_tp-PV_record_pl;
% r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname21,'/r2hill_record.txt'));
% r2hill_mean=mean(r2hill_record);
% hill21=sqrt(r2hill_mean);
% 
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname22,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname22,'/PV_record_tp.txt'));
% PV_rlt22=PV_record_tp-PV_record_pl;
% r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname22,'/r2hill_record.txt'));
% r2hill_mean=mean(r2hill_record);
% hill22=sqrt(r2hill_mean);

hillmax=3.5;
fontsize=15;

figure(1);
set(gcf,'Position',[400,100,700,700],'color','w');
%set(gca,'Position',[400,100,800,800],'color','w');

%subplot(2,2,1);
PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/PV_record_pl.txt'));
PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/PV_record_tp.txt'));
AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/AE_record_pl.txt'));
AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/AE_record_tp.txt'));

PV_rlt=PV_record_pl-PV_record_tp;
PV_rlt=PV_rlt(:,1:3);
Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));

PV_rlt_con=zeros(length(PV_rlt),3);
PV_rlt_con_norm=zeros(length(PV_rlt),1);
PV_tp_norm=zeros(length(PV_rlt),1);
for i=1:length(PV_rlt)
    P1=[cosd(AE_record_tp(i,4)) sind(AE_record_tp(i,4)) 0; ...,
    -sind(AE_record_tp(i,4)) cosd(AE_record_tp(i,4)) 0; ...,
    0 0 1];

    P2=[1 0 0;...,
    0 cosd(AE_record_tp(i,3)) sind(AE_record_tp(i,3)); ...,
    0 -sind(AE_record_tp(i,3)) cosd(AE_record_tp(i,3))];

    P3=[coswf(i) sinwf(i) 0;...,
    -sinwf(i) coswf(i) 0;...,
    0 0 1];
    PV_rlt_con(i,:)=(P3*P2*P1*PV_rlt(i,:)')';
    PV_rlt_con_norm(i)=norm(PV_rlt_con(i,:));
    PV_tp_norm(i)=norm(PV_record_tp(i,:));
end

r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);

N=50;
countx=zeros(N,1);
for i=1:N
    countx(i)=length(find(thetax(i) < theta & theta <= thetax(i+1)));
end

hill=sqrt(r2hill_mean);
xt=-pi:0.1:2*pi;
Rt=hillmax;
plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);hold all;axis square;
plot([-hillmax hillmax],[0 0],'r--');
plot([0 0],[-hillmax hillmax],'r--');
plot((PV_rlt_con(:,1).^2+PV_rlt_con(:,2).^2).^(1/2)/hill,PV_rlt_con(:,3)/hill,'k.');
%plot((PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2)/hill,PV_rlt(:,3)/hill,'r.');
xlim([-hillmax hillmax]);
ylim([-hillmax hillmax]);
xlabel('$$\Vert(xr,yr)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
ylabel('$$zr~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
title('1999CE119 & 2004UP10','fontsize',fontsize);
