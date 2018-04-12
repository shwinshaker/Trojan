clear;

Dir='ServerMount';
DDir='CEL4L5';
fname='1999CE119_1Gyr';
%fname='2001FU172_1Gyr';

tag='di';

CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);

tpelCE=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/AE_record_tp.txt']);
plelCE=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/AE_record_pl.txt']);
NepelCE=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/AE_record_Nep.txt']);

tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
Nepel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/Nepel.txt']);

timeCE=CE_record(:,1);
time=tpel(:,1);
phiCE=mod((tpelCE(:,4)+tpelCE(:,5)+tpelCE(:,6))-(NepelCE(:,4)+NepelCE(:,5)+NepelCE(:,6)),360);
phi=mod((tpel(:,5)+tpel(:,6)+tpel(:,7))-(Nepel(:,5)+Nepel(:,6)+Nepel(:,7)),360);

figure;

plot(time,phi,'k.');hold on;
plot(timeCE,phiCE,'r.');