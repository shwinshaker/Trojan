%% CE location distribution
clear;
ffname1='RealPlutinosNpl';
%ffname1='RanPlutinos';
ffname2='RealTrojans';
fname11='1999CE119_1Gyr_40pl';
%fname11='2001FU172_1Gyr_40pl';

%fname11='1999CE119_10k';
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
di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/di_record_inout.txt'));


PV_rlt=PV_record_pl-PV_record_tp;
PV_rlt=PV_rlt(:,1:3);
dis=(PV_rlt(:,1).^2+PV_rlt(:,2).^2+PV_rlt(:,3).^2).^(1/2);
Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
hnorm=((2*pi)^2/(365.25)^2*AE_record_tp(:,1).*(1-AE_record_tp(:,2).^2)).^(1/2);
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

miu=1.94148e-12;%(2*pi)^2/(365.25)^2;
absdi=abs(di_record_inout.*hnorm./Rnorm./coswf/miu);
maxdi=max(absdi);

r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname11,'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);

hill=sqrt(r2hill_mean);
xt=-pi:0.1:2*pi;
Rt=hillmax;
plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);hold all;axis square;
plot([-hillmax hillmax],[0 0],'r--');
plot([0 0],[-hillmax hillmax],'r--');
plotx=(PV_rlt_con(:,1).^2+PV_rlt_con(:,2).^2).^(1/2)/hill;
ploty=PV_rlt_con(:,3)/hill;
colorindex=abs(log(absdi/maxdi));
%colorindexNormal=(1+colorindex/10).^(-2);
colorindexNormal=(-erf(colorindex/3-2)+1)/2;
data=zeros(length(plotx),3);
data(:,1)=plotx;
data(:,2)=ploty;
data(:,3)=colorindexNormal;
datasort=sortrows(data,3);
color=[datasort(:,3) zeros(length(datasort),1) 1-datasort(:,3)];
% for i=1:length(PV_rlt_con)
%     plot(datasort(i,1),datasort(i,2),'.','color',[datasort(i,3) 0 1-datasort(i,3)],'markersize',12);
% end
scatter(datasort(:,1),datasort(:,2),20,color,'filled');
t=-pi/2:0.01:pi/2;
r=0.01:0.01:3.5-0.001;
[tt,rr]=meshgrid(t,r);
%zz=(-(rr.^2/hillmax^3)+1/hillmax+acos(rr./hillmax)./rr).*sin(tt);
zz=2*(1-(rr./hillmax).^2).^(1/2)./rr.*sin(tt);

[xx,yy]=pol2cart(tt,rr);
maxzz=max(max(zz));
ev=exp(exp(0:0.1:log(log(1+maxzz)+1))-1)-1;
ev=[ev -ev];
ev=unique(ev);
ev=sort(ev);
contour(xx,yy,zz,ev,'w--','linewidth',2);

hold off;
xlim([-hillmax hillmax]);
ylim([-hillmax hillmax]);
xlabel('$$\Vert(xr,yr)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
ylabel('$$zr~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
title('1999CE119 & 2004UP10','fontsize',fontsize);


% subplot(2,2,2);
% xt=-pi:0.1:2*pi;
% Rt=hillmax;
% plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);hold all;axis square;
% plot([-hillmax hillmax],[0 0],'r--');
% plot([0 0],[-hillmax hillmax],'r--');
% plot((PV_rlt12(:,1).^2+PV_rlt12(:,2).^2).^(1/2)/hill12,PV_rlt12(:,3)/hill12,'k.');
% xlim([-hillmax hillmax]);
% ylim([-hillmax hillmax]);
% xlabel('$$\Vert(xr,yr)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$$zr~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% title('2001FU172 & 2004UP10','fontsize',fontsize);
% hold off;
% subplot(2,2,3);
% xt=-pi:0.1:2*pi;
% Rt=hillmax;
% plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);hold all;axis square;
% plot([-hillmax hillmax],[0 0],'r--');
% plot([0 0],[-hillmax hillmax],'r--');
% plot((PV_rlt21(:,1).^2+PV_rlt21(:,2).^2).^(1/2)/hill21,PV_rlt21(:,3)/hill21,'k.');
% xlim([-hillmax hillmax]);
% ylim([-hillmax hillmax]);
% xlabel('$$\Vert(xr,yr)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$$zr~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% title('1999CE119 & 2006RJ103','fontsize',fontsize);
% hold off;
% subplot(2,2,4);
% xt=-pi:0.1:2*pi;
% Rt=hillmax;
% plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);hold all;axis square;
% plot([-hillmax hillmax],[0 0],'r--');
% plot([0 0],[-hillmax hillmax],'r--');
% plot((PV_rlt22(:,1).^2+PV_rlt22(:,2).^2).^(1/2)/hill22,PV_rlt22(:,3)/hill22,'k.');
% xlim([-hillmax hillmax]);
% ylim([-hillmax hillmax]);
% xlabel('$$\Vert(xr,yr)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$$zr~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
% title('2001FU172 & 2006RJ103','fontsize',fontsize);
% 
% hold off;