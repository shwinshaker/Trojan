clear;
ffname='RealPlutinos';
fname='1999CE119_1Gyr';

% tpel=load(['~/Documents/ServerMount/LAB Backup/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt']);
% plel=load(['~/Documents/ServerMount/LAB Backup/LAB/CE_realp/',ffname,'/',fname,'/plel.txt']);

tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt']);
plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/plel.txt']);

temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
if isempty(temp)
    ejectNo=size(tpel,1);
else
    ejectNo=temp;
end
clear temp;

tpel=tpel(1:ejectNo-1,:);
plel=plel(1:ejectNo-1,:);

aP=mean(plel(:,2));
eP=mean(plel(:,3));
IP=mean(plel(:,4))/180*pi;
aT=mean(tpel(:,2));
eT=mean(tpel(:,3));
IT=mean(tpel(:,4))/180*pi;

Ru=Fun_Ru(aP,aT,eP,eT,IP,IT);
disp(Ru)

Rk=Fun_Rk(aP,aT,1);

vp_mean=exp(log(1e-9):0.1:log(10));

figure;
fontsize=15;
loglog(vp_mean,(Ru/Rk)^2./vp_mean,'k-');
xlabel('$\overline{\nu_P}$','fontsize',fontsize,'interpreter','latex');
ylabel('$\widetilde{\nu_P}$','fontsize',fontsize,'interpreter','latex');

