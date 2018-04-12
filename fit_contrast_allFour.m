%% All Four
clear;

ffname='fit_record_contrast_201606';
fname='1999CE119_1Gyr';

name='1999CE119 & 2004UP10 1Gyr';

%tag='de';
tag='di';

load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));

load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_inout.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt'));

if strcmp(tag,'de')
    lim=0.001;
else
    lim=0.02;
end
xxlim=5.0;
yylim=10;
hill=sqrt(r2hill_record);
Row=2;
Column=2;

figure(1);
set(gcf,'Position',[400,100,600,600],'color','w');
annotation(gcf,'textbox','String',{[name,' ',tag]},'FontSize',12,'Position',[0.35 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))

subplot(Row,Column,1);
plot(CE_record(:,2)/hill,eval(strcat(tag,'_fit_perturb')),'k.');axis square;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(strcat(tag,'_fit_perturb')))))),'color','red');
grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' fit-perturb'));

subplot(Row,Column,2);
plot(CE_record(:,2)/hill,eval(strcat(tag,'_record_perturb')),'k.');axis square;
grid on;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(strcat(tag,'_record_perturb')))))),'color','red');
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' record-perturb'));

subplot(Row,Column,3);
plot(CE_record(:,2)/hill,eval(strcat(tag,'_fit_inout')),'k.');axis square;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(strcat(tag,'_fit_inout')))))),'color','red');
grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' fit-inout'));

subplot(Row,Column,4);
plot(CE_record(:,2)/hill,eval(strcat(tag,'_record_inout')),'k.');axis square;
grid on;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(strcat(tag,'_record_inout')))))),'color','red');
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' record-inout'));