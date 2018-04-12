clear;
%%%%% random fit inout & perturb

ffname='RanPlutinos';
fname='1999CE119_40pl';

%name='1999CE119 & 2004UP10 10plutoMass random';
name='1999CE119 ran';
tag='di';

load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));

%load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ran_CE_',tag,'.txt']);
load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ran_CE_',tag,'_ptb.txt']);

if strcmp(tag,'de')
    lim=0.001;
else
    lim=0.1;
end
xxlim=5.0;
yylim=100 ;
hill=sqrt(r2hill_record);
Row=2;
Column=2;

% PolyFit
% [P,H]=polyfit(log(eval(['ran_CE_',tag,'(:,1)'])/hill),log(abs((eval(['ran_CE_',tag,'(:,2)'])-eval(['ran_CE_',tag,'_ptb(:,2)']))./eval(['ran_CE_',tag,'_ptb(:,2)']))),1);  
% R=corrcoef(log(eval(['ran_CE_',tag,'(:,1)'])/hill),log(abs((eval(['ran_CE_',tag,'(:,2)'])-eval(['ran_CE_',tag,'_ptb(:,2)']))./eval(['ran_CE_',tag,'_ptb(:,2)']))));

% x2=0.1:0.01:xxlim;
% y2=exp(polyval(P,log(x2)));

figure(1);
set(gcf,'Position',[400,100,600,600],'color','w');
annotation(gcf,'textbox','String',{[name,' ',tag]},'FontSize',12,'Position',[0.35 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))


% subplot(2,2,1);
% plot(eval(['ran_CE_',tag,'(:,1)'])/hill,eval(['ran_CE_',tag,'(:,2)']),'k.');axis square;
% set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(['ran_CE_',tag,'(:,2)']))))),'color','red');
% set(text(1.5,5/6*lim,strcat('N=',num2str(length(eval(['ran_CE_',tag]))))),'color','red');
% 
% grid on;
% xlim([0 xxlim]);
% ylim([-lim lim]);
% xlabel('CE distance /Rhill');
% ylabel('\deltainc in this CE /deg');
% title(strcat(' fit inout'));

subplot(2,2,2);
plot(eval(['ran_CE_',tag,'_ptb(:,1)'])/hill,eval(['ran_CE_',tag,'_ptb(:,2)']),'k.');axis square;
set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval(['ran_CE_',tag,'_ptb(:,2)']))))),'color','red');
grid on;
xlim([0 xxlim]);
ylim([-lim lim]);
xlabel('CE distance /Rhill');
ylabel('\deltainc in this CE /deg');
title(strcat(' fit perturb'));

% subplot(2,2,3)
% loglog(eval(['ran_CE_',tag,'(:,1)'])/hill,abs((eval(['ran_CE_',tag,'(:,2)'])-eval(['ran_CE_',tag,'_ptb(:,2)']))./eval(['ran_CE_',tag,'_ptb(:,2)'])),'k.',...,
%     x2,y2,'r-');
% set(text(0.15,3/6*yylim,strcat('ln(y)=',num2str(P(1)),'ln(x)+',num2str(P(2)))),'fontsize',10,'color','red');
% set(text(0.15,2/6*yylim,strcat('R^2=',num2str(R(1,2)*R(2,1)))),'fontsize',10,'color','red');
% grid on;axis square;
% xlabel('CE distance /Rhill');
% ylabel('abs(inout-perturb)/perturb');
% title(strcat('fit:  abs(inout-perturb)/perturb'));
% xlim([0.1 xxlim]); 
% ylim([1/yylim yylim]);

% 
% plot(CE_record(:,2)/hill,di_fit_inout(:)./di_record_inout(:),'k.');
% grid on;
% xlabel('CE distance /Rhill');
% ylabel('fit/record');
% title(strcat('     fit-inout/record-inout'));
% xlim([0.1 xxlim]);
% ylim([-yylim yylim]);
% axis square;