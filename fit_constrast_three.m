%%CE_record;de_fit;de_record;r2hill
clear;
ffname='fit_record_contrast_201606';

fname='1999CE119_1Gyr';

tag='di';

name=['1999CE119 & 2004UP10 1Gyr ',tag];

load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_inout.txt'));

if strcmp(tag,'de')
    %lim=0.001;
    lim=0.0002;
else
    lim=0.02;
    %lim=0.002;
end
xxlim=5.0;
yylim=10 ;
hill=sqrt(r2hill_record);
Row=2;
Column=3;



% PolyFit
% [P,H]=polyfit(log(CE_record(:,2)/hill),log(abs((de_fit_perturb(:)-de_record_perturb(:))./de_record_perturb(:))),1);  
% R=corrcoef(log(CE_record(:,2)/hill),log(abs((de_fit_perturb(:)-de_record_perturb(:))./de_record_perturb(:))));

figure(1);
set(gcf,'Position',[400,100,900,600],'color','w');
annotation(gcf,'textbox','String',{name},'FontSize',12,'Position',[0.40 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))

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
plot(CE_record(:,2)/hill,eval(strcat(tag,'_record_inout')),'k.');axis square;
grid on;
set(text(1.5,5/6*lim,strcat('N=',num2str(length(CE_record)))),'color','red');
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(strcat(tag,'_record_inout')))))),'color','red');
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' record-inout'));

%saveas(1,name,'bmp');
% print(gcf,'-dbitmap',name);

% subplot(Row,Column,4);
% defitpoly=exp(polyval(P,log(CE_record(:,2)/hill)));
% %figure(5);
% %loglog(CE_record(:,2)/hill,defitpoly(:),'k.');
% plot(CE_record(:,2)/hill,de_fit_perturb(:)./(1-defitpoly(:)),'k.');axis square;
% grid on;
% set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(de_fit_perturb(:)./(1-defitpoly(:)))))),'color','red');
% xlim([0 xxlim]);
% ylim([-lim lim]);
% xlabel('CE distance /Rhill');
% ylabel('\deltainc in this CE /deg');
% title(strcat(' Modefied fit '));
% 
% subplot(Row,Column,3);
% x2=0.1:0.01:xxlim;
% y2=exp(polyval(P,log(x2)));
% loglog(CE_record(:,2)/hill,abs((de_fit_perturb(:)-de_record_perturb(:))./de_record_perturb(:)),'k.',...,
%     x2,y2,'r-');
% set(text(0.15,3/4*yylim,strcat('ln(y)=',num2str(P(1)),'ln(x)+',num2str(P(2)))),'fontsize',15,'color','red');
% set(text(0.15,2/4*yylim,strcat('R^2=',num2str(R(1,2)*R(2,1)))),'fontsize',15,'color','red');
% grid on;axis square;
% %loglog(CE_record(:,2)/hill,abs(de_fit_perturb(:)./de_record_perturb(:)),'ko');
% xlabel('CE distance /Rhill');
% ylabel('(fit-record)/record');
% title(strcat('      (fit-record)/record'));
% xlim([0.1 xxlim]); 
% ylim([1/yylim yylim]);
% % print(gcf,'-dbitmap');

subplot(Row,Column,4);
%plot(CE_record(:,2)/hill,(de_fit_perturb(:)-de_record_perturb(:))./de_record_perturb(:),'ko');
plot(CE_record(:,2)/hill,eval(strcat(tag,'_fit_perturb'))./eval(strcat(tag,'_record_perturb')),'k.');
grid on;
xlabel('CE distance /Rhill');
ylabel('fit/record');
title(strcat('     fit-perturb/record-perturb'));
xlim([0.1 xxlim]);
ylim([-yylim yylim]);
axis square;
% print(gcf,'-dbitmap');

%figure(2);
subplot(Row,Column,5);
% set (gcf,'Position',[400,100,400,400], 'color','w')
plot(CE_record(:,2)/hill,eval(strcat(tag,'_record_perturb'))./eval(strcat(tag,'_record_inout')),'k.');
grid on;
xlabel('CE distance /Rhill');
ylabel('perturb/inout');
title(strcat('     record-perturb/record-inout'));
xlim([0.1 xxlim]);
ylim([-yylim yylim]);
axis square;


% % The deviation to 1 of fit/record
% plot(CE_record(:,2)/hill,(de_fit_perturb(:)./de_record_perturb(:)-1).^2,'k.');
% grid on;
% xlabel('CE distance /Rhill');
% ylabel('Deviation');
% title(strcat('     (fit-perturb/record-perturb-1)^2'));
% xlim([0.1 xxlim]);
% ylim([-yylim yylim]);
% axis square;