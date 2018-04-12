%% total i change
clear;

ffname='fit_record_contrast_201606';
fname='1999CE119_forgedMass';

%tag=2;

di_perturb_smlstp=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/di_perturb_smlstp.txt'));
di_inout_smlstp=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/di_inout_smlstp.txt'));
dis_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/dis_record.txt'));
r2hill_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
hill=sqrt(r2hill_record);
xxlim=5.0;
yylim1=5e-2;
yylim=1000;
N=length(dis_record)-1;


 
% PolyFit
[P,H]=polyfit(log(dis_record/hill),log(abs((di_inout_smlstp-di_perturb_smlstp)./di_perturb_smlstp)),1);  
R=corrcoef(log(dis_record/hill),log(abs((di_inout_smlstp-di_perturb_smlstp)./di_perturb_smlstp)));

x2=0.1:0.01:xxlim;
y2=exp(polyval(P,log(x2)));

% figure(1);
% set(gcf,'Position',[400,100,1000,500],'color','w');
% for i=1:N
%     
% subplot(1,2,1);
% 
% plot(dis_record(i)/hill,di_perturb_smlstp(i),'b.');hold on;
% xlim([0 xxlim]);ylim([-5.5e-5 5.5e-5]);
% axis square;grid on;
% plot(dis_record(i)/hill,di_inout_smlstp(i),'r.');hold on;
% pause(0.05);
% 
% subplot(1,2,2);
% 
% clr=i/N;
% loglog(dis_record(i)/hill,-(di_inout_smlstp(i)-di_perturb_smlstp(i))./di_perturb_smlstp(i),'.','color',[clr 0 1-clr]);
% hold on;
% xlim([0.1 xxlim]); 
% ylim([1/yylim yylim]);
% grid on;axis square;
% pause(0.05);
% end

figure(2);
set(gcf,'Position',[400,100,1000,500],'color','w');

subplot(1,2,1);
plot(dis_record/hill,di_perturb_smlstp,'b.');hold on;

xlim([0 xxlim]);ylim([-yylim1 yylim1]);
axis square;grid on;
plot(dis_record/hill,di_inout_smlstp,'r.');
legend('Perturb','Inout');
xlabel('CE distance /Rhill');
ylabel('\deltainc in the small step /deg');
title(strcat(' Contrast'));
text(0.6*xxlim,0.7*yylim1,['Tot \Deltainc: ',num2str(sum(di_inout_smlstp),'%.2e'),' Deg'],'color','r');
text(0.6*xxlim,0.6*yylim1,['Tot \Deltainc: ',num2str(sum(di_perturb_smlstp),'%.2e'),' Deg'],'color','b');

subplot(1,2,2);

loglog(dis_record/hill,abs((di_inout_smlstp-di_perturb_smlstp)./di_perturb_smlstp),'k.');
hold on;
xlim([0.1 xxlim]); 
ylim([1/yylim yylim]);
grid on;axis square;

loglog(x2,y2,'r-');
set(text(0.15,3/4*yylim,strcat('ln(y)=',num2str(P(1)),'ln(x)',num2str(P(2)))),'fontsize',15,'color','red');
set(text(0.15,2/4*yylim,strcat('R^2=',num2str(R(1,2)*R(2,1)))),'fontsize',15,'color','red');
xlabel('CE distance /Rhill');
ylabel('abs(Inout-Perturb)/Perturb');
title(strcat('      abs(Inout-Perturb)/Perturb'));



% loglog(dis_record(1:N)/hill,(di_inout_smlstp(1:N)-di_perturb_smlstp(1:N))./di_perturb_smlstp(1:N),'k.');
