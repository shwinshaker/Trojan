%%MassTest
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')
clear;
Dir='ServerMount'; 
fontsize=15;
PluName='1999CE119';
di_name='de_record_inout';
%di_name='de_record_inout';
RatioList1=1./exp(log(1):0.1:log(1000))';
RatioList1=sort(RatioList1);
RatioList2=exp(log(1):0.1:log(1000))';
RatioList=[RatioList1;RatioList2];

fname=cell(length(RatioList),1);
for i=1:length(RatioList1)
    Ratio=RatioList1(i);
    fname{i}=[PluName,'_',num2str(sprintf('%.4f',Ratio)),'MP'];
end
for i=1:length(RatioList2)
    Ratio=RatioList2(i);
    fname{length(RatioList1)+i}=[PluName,'_',num2str(sprintf('%.2f',Ratio)),'MP'];
end

%%%%%%%%%%%%% Revise 1&2 %%%%%%%%%%%%%
barrier=5.0;
iF=find(RatioList2<=barrier,1,'last');
RatioList1=[RatioList1;RatioList2(1:iF)];
RatioList2(1:iF)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%
TimesList=zeros(length(RatioList),1);
SumList=zeros(length(RatioList),1);
StdList=zeros(length(RatioList),1);
RList=zeros(length(RatioList),1);
EndTime=zeros(length(RatioList),1);
for i=1:length(RatioList)
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/MassTest2/',fname{i},'/',di_name,'.txt']);
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/MassTest2/',fname{i},'/CE_record.txt']);
    EndTime(i)=CE_record(end,1)/365;
    CE_times=length(di_record_inout);
    absSumdi=abs(sum(di_record_inout));
    RList(i)=sum((di_record_inout-mean(di_record_inout)).^2);
    Stddi=(var(di_record_inout))^(1/2);
    %Stddi=(RList(i)/CE_times)^(1/2);
    TimesList(i)=CE_times;
    SumList(i)=absSumdi;
    StdList(i)=Stddi;
end

[PT1,HT1]=polyfit(log(RatioList1),log(TimesList(1:length(RatioList1))),1);  
RT1=corrcoef(log(RatioList1),log(TimesList(1:length(RatioList1))));
Timesfit1=exp(polyval(PT1,log(RatioList1)));

[PT2,HT2]=polyfit(log(RatioList2),log(TimesList(length(RatioList1)+1:end)),1);  
RT2=corrcoef(log(RatioList2),log(TimesList(length(RatioList1)+1:end)));
Timesfit2=exp(polyval(PT2,log(RatioList2)));

[PSum1,HSum1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1))),1);  
RSum1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1))));
Sumfit1=exp(polyval(PSum1,log(RatioList1)));

[PSum2,HSum2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end)),1);  
RSum2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end)));
Sumfit2=exp(polyval(PSum2,log(RatioList2)));

[PStd1,HStd1]=polyfit(log(RatioList1),log(StdList(1:length(RatioList1))),1);  
RStd1=corrcoef(log(RatioList1),log(StdList(1:length(RatioList1))));
Stdfit1=exp(polyval(PStd1,log(RatioList1)));

[PStd2,HStd2]=polyfit(log(RatioList2),log(StdList(length(RatioList1)+1:end)),1);  
RStd2=corrcoef(log(RatioList2),log(StdList(length(RatioList1)+1:end)));
Stdfit2=exp(polyval(PStd2,log(RatioList2)));

[PR1,HR1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1))),1);  
RR1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1))));
Rfit1=exp(polyval(PR1,log(RatioList1)));

[PR2,HR2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end)),1);  
RR2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end)));
Rfit2=exp(polyval(PR2,log(RatioList2)));

end


linewidth=2;
xxlim=[1e-3 1e3];
xxtick=power(10,-3:1:3);
xp=1.5e-3;

figure;

set(gcf,'Position',[400,100,600,700],'color','w');
row=2;
col=2;
subplot(row,col,1);
loglog(RatioList,TimesList,'k.');hold all;
loglog(RatioList1,Timesfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Timesfit2,'r-','linewidth',linewidth);
Timesfit1_extend=exp(polyval(PT1,log(RatioList2)));
loglog(RatioList2,Timesfit1_extend,'r-.','linewidth',linewidth);
yylim=[10 1e6];
ylim(yylim);
xlim(xxlim);
plot([barrier barrier],yylim,'r--');
hold off;

xlabel('$M/M_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'xTick',xxtick);

set(text(xp,5/10*yylim(2),['$$\ln{y} = ',num2str(PT1(1),'%.4f'),'\,\ln{x}+',num2str(PT1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,2/10*yylim(2),['$$R^2 = ',num2str(RT1(1,2)*RT1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(1e-2,5/1000*yylim(2),['$$\ln{y} = ',num2str(PT2(1),'%.4f'),'\,\ln{x}+',num2str(PT2(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');
% set(text(1e-2,3/1000*yylim(2),['$$R^2 = ',num2str(RT2(1,2)*RT2(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');
title('(A)');

subplot(row,col,2);
loglog(RatioList,SumList,'k.');hold all;
loglog(RatioList1,Sumfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Sumfit2,'r-','linewidth',linewidth);
Sumfit1_extend=exp(polyval(PSum1,log(RatioList2)));
loglog(RatioList2,Sumfit1_extend,'r-.','linewidth',linewidth);

xlabel('$M/M_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$A_I=\sum \Delta I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
yylim=[1e-7 1e3];
ylim(yylim);
xlim(xxlim);

plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSum1(1),'%.4f'),'\,\ln{x}',num2str(PSum1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.5/10*yylim(2),['$$R^2 = ',num2str(RSum1(1,2)*RSum1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(1e-2,5/1e7*yylim(2),['$$\ln{y} = ',num2str(PSum2(1),'%.4f'),'\,\ln{x}',num2str(PSum2(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');
% set(text(1e-2,10/1e8*yylim(2),['$$R^2 = ',num2str(RSum2(1,2)*RSum2(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');
title('(B)');


% subplot(2,1,2);
% loglog(RatioList,StdList,'k.-');hold on;
% loglog(RatioList,Stdfit,'r.-');
subplot(row,col,3);
loglog(RatioList,RList,'k.');hold all;
loglog(RatioList1,Rfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Rfit2,'r-','linewidth',linewidth);
Rfit1_extend=exp(polyval(PR1,log(RatioList2)));
loglog(RatioList2,Rfit1_extend,'r-.','linewidth',linewidth);

yylim=[1e-10 1e5];
ylim(yylim);
xlim(xxlim);

%yylim=get(gca,'ylim');
plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PR1(1),'%.4f'),'\,\ln{x}',num2str(PR1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,1/100*yylim(2),['$$R^2 = ',num2str(RR1(1,2)*RR1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(1e-1,1/1e10*yylim(2),['$$\ln{y} = ',num2str(PR2(1),'%.4f'),'\,\ln{x}',num2str(PR2(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');
% set(text(1e-1,2/1e11*yylim(2),['$$R^2 = ',num2str(RR2(1,2)*RR2(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','blue');

xlabel('$M/M_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
title('(C)');


subplot(row,col,4);
loglog(RatioList,EndTime,'k.');hold on;
yylim=[1e5 1e10];
ylim(yylim);
xlim(xxlim);

xxlim=get(gca,'xlim');
plot([barrier barrier],yylim,'r--');
%plot(xxlim,[1e9 1e9],'r--');
xlabel('$M/M_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Evolution~Time~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
title('(D)');
set(gca,'xTick',xxtick);

% xlabel('$M/M_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$M_0$','fontsize',fontsize,'Interpreter','latex');
% yylim=get(gca,'ylim');
% set(text(2e-3,2/10*yylim(2),['$$\ln{y} = ',num2str(P4(1),'%.4f'),'\,\ln{x}+',num2str(P4(2),'%.4f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');


