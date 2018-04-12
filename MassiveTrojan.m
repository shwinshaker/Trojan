clear;

Dir='ServerMount';
DDir='MassiveTrojan2';

fname='1999CE119_1Gyr_40pl';

tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel43.txt']);
absDel=max(abs(tpel(:,3)-mean(tpel(:,3))));

