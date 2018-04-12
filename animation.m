%% animate
clear;

ffname='fit_record_contrast_201606';
fname='1999CE119_1Gyr';

name='1999CE119 & 2004UP10 1Gyr';

tag='di';
load(strcat('~/Documents/swiftdata/Trojan/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
load(strcat('~/Documents/swiftdata/Trojan/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
load(strcat('~/Documents/swiftdata/Trojan/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt'));


lim=0.02;
xxlim=5.0;
yylim=100;
hill=sqrt(r2hill_record);
Row=2;
Column=2;

figure(1);
digits(4);

maxtime=max(CE_record(:,1));

for i=1:length(CE_record)

clr=CE_record(i,1)/maxtime;
plot(CE_record(i,2)/hill,eval(strcat(tag,'_record_inout(i)')),'.','color',[clr 0 1-clr]);hold on;
axis square;
grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(name);
h=text(1.5,3/4*lim,strcat('sum=',char(vpa(sum(eval(strcat(tag,'_record_inout(1:i)')))))),'color','red','fontsize',20);
pause(0.01);
if i~=length(CE_record) 
    set(h,'visible','off');
end
end