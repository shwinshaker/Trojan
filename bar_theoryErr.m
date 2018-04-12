%% CE location distribution
clear;
% ffname='RealPlutinosNpl';
% %ffname='RanPlutinos';
% %ffname='RealPlutinos';
% %ffname2='RealTrojans';
% fname='1999CE119_1Gyr_40pl';
% %fname='1999CE119&2006RJ103_1Gyr_40pl';
% %fname='2001FU172_1Gyr_40pl';
% %fname='2001FU172&2006RJ103_1Gyr_40pl';
% %fname='2000FB8_1Gyr';
% %fname='1998HK151_1Gyr';
% %fname='1999CE119_10k';
ffname='RealPlutinosNpl';
fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};


tag='';
%tag='ran';

if strcmp(tag,'ran')
    difile='di_fit_perturb.txt';
    PVtpfile='ran_PV_record_tp.txt';
    PVplfile='ran_PV_record_pl.txt';
    AEtpfile='ran_AE_record_tp.txt';
    AEplfile='ran_AE_record_pl.txt';    
else
    difile='di_record_perturb.txt';
    PVtpfile='PV_record_tp.txt';
    PVplfile='PV_record_pl.txt';
    AEtpfile='AE_record_tp.txt';
    AEplfile='AE_record_pl.txt'; 
end

hillmax=3.5;
fontsize=15;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);
hill=sqrt(r2hill_mean);

PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',PVplfile));
PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',PVtpfile));
AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',AEplfile));
AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',AEtpfile));
di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',difile));

if strcmp(tag,'ran')
    ierror=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt'));
    PV_record_pl(ierror,:)=[];
    PV_record_tp(ierror,:)=[];
    AE_record_pl(ierror,:)=[];
    AE_record_tp(ierror,:)=[];
end
PV_rlt=PV_record_pl-PV_record_tp;
P_rlt=PV_rlt(:,1:3);
V_rlt=PV_rlt(:,4:6);
dishill=(P_rlt(:,1).^2+P_rlt(:,2).^2+P_rlt(:,3).^2).^(1/2)/hill;
Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
Rpnorm=(PV_record_pl(:,1).^2+PV_record_pl(:,2).^2+PV_record_pl(:,3).^2).^(1/2);

Vrnorm=(V_rlt(:,1).^2+V_rlt(:,2).^2+V_rlt(:,3).^2).^(1/2);
hnorm=((2*pi)^2/(365.25)^2*AE_record_tp(:,1).*(1-AE_record_tp(:,2).^2)).^(1/2);
sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));

PV_rlt_con=zeros(length(P_rlt),3);
PV_rlt_con_norm=zeros(length(P_rlt),1);
P_tp_norm=zeros(length(P_rlt),1);

for i=1:length(P_rlt)
    P1=[cosd(AE_record_tp(i,4)) sind(AE_record_tp(i,4)) 0; ...,
    -sind(AE_record_tp(i,4)) cosd(AE_record_tp(i,4)) 0; ...,
    0 0 1];

    P2=[1 0 0;...,
    0 cosd(AE_record_tp(i,3)) sind(AE_record_tp(i,3)); ...,
    0 -sind(AE_record_tp(i,3)) cosd(AE_record_tp(i,3))];

    P3=[coswf(i) sinwf(i) 0;...,
    -sinwf(i) coswf(i) 0;...,
    0 0 1];
    PV_rlt_con(i,:)=(P3*P2*P1*P_rlt(i,:)')';
    PV_rlt_con_norm(i)=norm(PV_rlt_con(i,:));
    P_tp_norm(i)=norm(PV_record_tp(i,:));
end

theta=atan(PV_rlt_con(:,3)./(PV_rlt_con(:,1).^2+PV_rlt_con(:,2).^2).^(1/2));

miu=1.94148e-12;
di_nor=di_record_inout/180*pi.*hnorm./Rnorm/miu;
dis=dishill*hill;
maxdi=max(di_nor);

C=mean(abs(2*miu*Rnorm./hnorm));
disp(fname{isub});
disp(C);
disp([mean(Rpnorm),mean(Rnorm),mean(hnorm),mean(abs(coswf)),mean(Vrnorm)]);

%%pick the dis in range
hillll=0.0;
hillul=3.5;

di_norx=di_nor(hillll< dishill & dishill <hillul);

disx=dis(hillll< dishill & dishill <hillul);
thetax=theta(hillll< dishill & dishill <hillul);
Vrnormx=Vrnorm(hillll< dishill & dishill <hillul);
coswfx=coswf(hillll< dishill & dishill <hillul);

%plot(thetax,di_norx,'k.');hold all;
%plot(disx/hill,di_norx,'k.');hold all;
%dis0=2.55*hill;
%t=-pi/2:0.01:pi/2;

dismax=hillmax*hill;
%zz=(-(dis0^2/dismax^3)+1/dismax+acos(dis0/dismax)/dis0)*sin(t);
%zz=(-(disx.^2/dismax^3)+1/dismax+acos(disx./dismax)./disx).*sin(thetax);
zz=2*(1-(disx/dismax).^2).^(1/2)./disx.*sin(thetax)./Vrnormx.*coswfx;
%plot(thetax,zz,'r.');
%plot(disx/hill,zz,'r.');

rlterr=abs((zz-di_norx)./di_norx);
% rlterrabs=abs(rlterr);
% rlterrp=rlterr(rlterr>=0);
% thetaxp=thetax(rlterr>=0);
% rlterrn=rlterr(rlterr<0);
% thetaxn=thetax(rlterr<0);

rlterrsort=sort(rlterr);
mark=rlterrsort(floor(length(rlterrsort)*0.90));
% 
% subplot(1,2,1);
% semilogy(thetax/pi*180,rlterr,'k.');hold all;%axis square;
% 
% plot([-90 90],[median(rlterr) median(rlterr)],'r-');
% plot([-90 90],[mark mark],'r-','linewidth',3);
% hold off;
% xlim([-90 90]);
% set(gca,'xTick',-90:30:90);
% xlabel('$\theta~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$Error$','fontsize',fontsize,'Interpreter','latex');

subplot(2,2,isub);
semilogy(disx/hill,rlterr,'k.');hold all;%axis square;

plot([-90 90],[median(rlterr) median(rlterr)],'r-');
plot([-90 90],[mark mark],'r-','linewidth',3);
hold off;
xlim([0 3.5]);
ylim([1e-8 1e4]);
title(titlename{isub},'fontsize',fontsize);

xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Err$','fontsize',fontsize,'Interpreter','latex');


% [xx,yy]=pol2cart(tt,rr);
% maxzz=max(max(zz));
% ev=exp(exp(0:0.1:log(log(1+maxzz)+1))-1)-1;
% ev=[ev -ev];
% ev=unique(ev);
% ev=sort(ev);
% contour(xx,yy,zz,ev,'w--','linewidth',2);
hold off;

end
% xlim([-hillmax hillmax]);
% ylim([-hillmax hillmax]);

% title('1999CE119 & 2004UP10','fontsize',fontsize);