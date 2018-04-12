clear;

ffname='fit_record_contrast_201606';
fname='1999CE119_1Gyr';

tag='de';
name='1999CE119 & 2004UP10 1Gyr';


load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));

load(['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt']);
load(['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_inout.txt']);

if strcmp(tag,'de')
    lim=0.001;
else
    lim=0.02;
end
xxlim=5.0;
yylim=100 ;
hill=sqrt(r2hill_record);
Row=2;
Column=2;

% PolyFit
[P,H]=polyfit(log(CE_record(:,2)/hill),log(abs((eval([tag,'_fit_inout'])-eval([tag,'_fit_perturb']))./eval([tag,'_fit_perturb']))),1);  
R=corrcoef(log(CE_record(:,2)/hill),log(abs((eval([tag,'_fit_inout'])-eval([tag,'_fit_perturb']))./eval([tag,'_fit_perturb']))));

x2=0.1:0.01:xxlim;
y2=exp(polyval(P,log(x2)));

figure(1);
set(gcf,'Position',[400,100,600,600],'color','w');
annotation(gcf,'textbox','String',{[name,' ',tag]},'FontSize',12,'Position',[0.35 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))

subplot(2,2,1);
plot(CE_record(:,2)/hill,eval([tag,'_fit_inout']),'k.');axis square;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval([tag,'_fit_inout']))))),'color','red');
set(text(1.5,5/6*lim,strcat('N=',num2str(length(CE_record)))),'color','red');

grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' fit inout'));

subplot(2,2,2);
plot(CE_record(:,2)/hill,eval([tag,'_fit_perturb']),'k.');axis square;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval([tag,'_fit_perturb']))))),'color','red');
grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' fit perturb'));

subplot(2,2,3)
loglog(CE_record(:,2)/hill,abs((eval([tag,'_fit_inout'])-eval([tag,'_fit_perturb']))./eval([tag,'_fit_perturb'])),'k.',...,
    x2,y2,'r-');
set(text(0.15,3/6*yylim,strcat('ln(y)=',num2str(P(1)),'ln(x)+',num2str(P(2)))),'fontsize',10,'color','red');
set(text(0.15,2/6*yylim,strcat('R^2=',num2str(R(1,2)*R(2,1)))),'fontsize',10,'color','red');
grid on;axis square;
xlabel('CE distance /Rhill');
ylabel('abs(inout-perturb)/perturb');
title(strcat('fit:  abs(inout-perturb)/perturb'));
xlim([0.1 xxlim]); 
ylim([1/yylim yylim]);



% 
% plot(CE_record(:,2)/hill,di_fit_inout(:)./di_record_inout(:),'k.');
% grid on;
% xlabel('CE distance /Rhill');
% ylabel('fit/record');
% title(strcat('     fit-inout/record-inout'));
% xlim([0.1 xxlim]);
% ylim([-yylim yylim]);
% axis square;