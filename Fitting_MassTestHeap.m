%%MassTest
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')
clear;
Dir='ServerMount'; 

% PluName='1999CE119';
% DDir='MassTest2';
% DDir='MassTest3';
% PluName='2001FU172';
% PluName='1999CE119_2006RJ103';
% DDir='MassTest4';
PluName='2001FU172_2006RJ103';
DDir='MassTest5';

% PluName='1999CE119';
% DDir='MassTest2_denser';

di_name='di_record_inout';
de_name='de_record_inout';
da_name='da_record_inout';

RatioList1=1./exp(log(1):0.1:log(1000))';
RatioList1=sort(RatioList1);
RatioList2=exp(log(1):0.1:log(1000))';
RatioList=[RatioList1;RatioList2];

% RatioList1=exp(log(1e-3):0.1:log(1))';
% %RatioList1=sort(RatioList1);
% RatioList2=exp(log(1):0.1:log(1e3))';
% RatioList=[RatioList1;RatioList2];

% RatioList1=exp(log(1e-3):0.05:log(1))';
% %RatioList1=sort(RatioList1);
% RatioList2=exp(log(1):0.05:log(1e3))';
% RatioList=[RatioList1;RatioList2];

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
if strcmp(PluName,'1999CE119') 
    barrier=14.8694;
elseif strcmp(PluName,'1999CE119_2006RJ103') 
    barrier=10.0;
else
    barrier=30.0;
end

iF=find(RatioList2<=barrier,1,'last');
RatioList1=[RatioList1;RatioList2(1:iF)];
RatioList2(1:iF)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%
MinDisList=zeros(length(RatioList),1);

aPL=zeros(length(RatioList),1);
IPL=zeros(length(RatioList),1);
ePL=zeros(length(RatioList),1);
aTL=zeros(length(RatioList),1);
ITL=zeros(length(RatioList),1);

TimesList=zeros(length(RatioList),1);

MaxList=zeros(length(RatioList),3);
MaxTime=zeros(length(RatioList),3);
SumList=zeros(length(RatioList),3);
StdList=zeros(length(RatioList),3);
RList=zeros(length(RatioList),3);
EjectTime=zeros(length(RatioList),1);
absDelList=zeros(length(RatioList),3);
MaxCumsumList=zeros(length(RatioList),3);
for i=1:length(RatioList)
    disp(fname{i});
    
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/plel.txt']);
    
    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if isempty(temp)
        ejectNo=size(tpel,1);
    else
        ejectNo=temp;
    end
    clear temp;
    ejecttime=tpel(ejectNo,1);
    
    tpel=tpel(1:ejectNo-1,:);
    plel=plel(1:ejectNo-1,:);
    
    aPL(i)=mean(plel(:,2));
    IPL(i)=mean(plel(:,4))/180*pi;
    ePL(i)=mean(plel(:,3));
    aTL(i)=mean(tpel(:,2));
    ITL(i)=mean(tpel(:,4))/180*pi;
    
    %absDelList(i,1)=abs(tpel(end,4)-tpel(1,4));
    absDelList(i,1)=max(abs(tpel(:,4)-mean(tpel(:,4))))/180*pi;
    %absDelList(i,2)=abs(tpel(end,3)-tpel(1,3));
    absDelList(i,2)=max(abs(tpel(:,3)-mean(tpel(:,3))));
    absDelList(i,3)=max(abs(tpel(:,2)-mean(tpel(:,2))));
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);
%     da_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',da_name,'.txt']);

    CEtime=CE_record(1:ejectCENo,1)/365.25;
    di_record_inout=di_record_inout(1:ejectCENo)/180*pi;
    de_record_inout=de_record_inout(1:ejectCENo);
%     da_record_inout=da_record_inout(1:ejectCENo);

    di_sum=cumsum(di_record_inout);
    de_sum=cumsum(de_record_inout);
%     da_sum=cumsum(da_record_inout);
    
    EjectTime(i)=ejecttime/365.25;
    TimesList(i)=length(di_record_inout);
    
    MinDisList(i)=min(CE_record(:,2));
    
    [MaxList(i,1),ind]=max(abs(di_record_inout));
    %MaxTime(i,1)=CEtime(abs(di_record_inout)==max(abs(di_record_inout)));
    MaxTime(i,1)=CEtime(ind);

    SumList(i,1)=abs(sum(di_record_inout));
    RList(i,1)=sum((di_record_inout-mean(di_record_inout)).^2);
    StdList(i,1)=std(di_record_inout);
    MaxCumsumList(i,1)=max(di_sum)-min(di_sum);
    
    [MaxList(i,2),ind]=max(abs(de_record_inout));
    MaxTime(i,2)=CEtime(ind);
    
    SumList(i,2)=abs(sum(de_record_inout));
    RList(i,2)=sum((de_record_inout-mean(de_record_inout)).^2);
    StdList(i,2)=std(de_record_inout);
    MaxCumsumList(i,2)=max(de_sum)-min(de_sum);
    
%     [MaxList(i,3),ind]=max(abs(da_record_inout));
%     MaxTime(i,3)=CEtime(ind);
%     
%     SumList(i,3)=abs(sum(da_record_inout));
%     RList(i,3)=sum((da_record_inout-mean(da_record_inout)).^2);
%     StdList(i,3)=std(da_record_inout);
%     MaxCumsumList(i,3)=max(da_sum)-min(da_sum);
end

aP=mean(aPL);
aT=mean(aTL);
eP=mean(ePL);
IP=mean(IPL);
IT=mean(ITL);

[PT1,HT1]=polyfit(log(RatioList1),log(TimesList(1:length(RatioList1),1)),1);  
RT1=corrcoef(log(RatioList1),log(TimesList(1:length(RatioList1),1)));
Timesfit1=exp(polyval(PT1,log(RatioList1)));

[PT2,HT2]=polyfit(log(RatioList2),log(TimesList(length(RatioList1)+1:end,1)),1);  
RT2=corrcoef(log(RatioList2),log(TimesList(length(RatioList1)+1:end,1)));
Timesfit2=exp(polyval(PT2,log(RatioList2)));

[PSum1,HSum1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),1)),1);  
RSum1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),1)));
Sumfit1=exp(polyval(PSum1,log(RatioList1)));

[PSum2,HSum2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,1)),1);  
RSum2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,1)));
Sumfit2=exp(polyval(PSum2,log(RatioList2)));

%%% Cum Sum Max
[PCumsum1,HCumsum1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)),1);  
RCumsum1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)));
Cumsumfit1=exp(polyval(PCumsum1,log(RatioList1)));

[PCumsum2,HCumsum2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)),1);  
RCumsum2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)));
Cumsumfit2=exp(polyval(PCumsum2,log(RatioList2)));

%%%
[PStd1,HStd1]=polyfit(log(RatioList1),log(StdList(1:length(RatioList1),1)),1);  
RStd1=corrcoef(log(RatioList1),log(StdList(1:length(RatioList1),1)));
Stdfit1=exp(polyval(PStd1,log(RatioList1)));

[PStd2,HStd2]=polyfit(log(RatioList2),log(StdList(length(RatioList1)+1:end,1)),1);  
RStd2=corrcoef(log(RatioList2),log(StdList(length(RatioList1)+1:end,1)));
Stdfit2=exp(polyval(PStd2,log(RatioList2)));

[PR1,HR1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),1)),1);  
RR1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),1)));
Rfit1=exp(polyval(PR1,log(RatioList1)));

[PR2,HR2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,1)),1);  
RR2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,1)));
Rfit2=exp(polyval(PR2,log(RatioList2)));

%%de
[PSumde1,HSumde1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),2)),1);  
RSumde1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),2)));
Sumfitde1=exp(polyval(PSumde1,log(RatioList1)));

[PSumde2,HSumde2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,2)),1);  
RSumde2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,2)));
Sumfitde2=exp(polyval(PSumde2,log(RatioList2)));

%% Cum Sum Max
[PCumsumde1,HCumsumde1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)),1);  
RCumsumde1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)));
Cumsumfitde1=exp(polyval(PCumsumde1,log(RatioList1)));

[PCumsumde2,HCumsumde2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)),1);  
RCumsumde2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)));
Cumsumfitde2=exp(polyval(PCumsumde2,log(RatioList2)));


[PRde1,HRde1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),2)),1);  
RRde1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),2)));
Rfitde1=exp(polyval(PRde1,log(RatioList1)));

[PRde2,HRde2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,2)),1);  
RRde2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,2)));
Rfitde2=exp(polyval(PRde2,log(RatioList2)));

%%da
[PSumda1,HSumda1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),3)),1);  
RSumda1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),3)));
Sumfitda1=exp(polyval(PSumda1,log(RatioList1)));

[PSumda2,HSumda2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,3)),1);  
RSumda2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,3)));
Sumfitda2=exp(polyval(PSumda2,log(RatioList2)));

%% Cum Sum Max
[PCumsumda1,HCumsumda1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),3)),1);  
RCumsumda1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),3)));
Cumsumfitda1=exp(polyval(PCumsumda1,log(RatioList1)));

[PCumsumda2,HCumsumda2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,3)),1);  
RCumsumda2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,3)));
Cumsumfitda2=exp(polyval(PCumsumda2,log(RatioList2)));


[PRda1,HRda1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),3)),1);  
RRda1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),3)));
Rfitda1=exp(polyval(PRda1,log(RatioList1)));

[PRda2,HRda2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,3)),1);  
RRda2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,3)));
Rfitda2=exp(polyval(PRda2,log(RatioList2)));

end


linewidth=2;
fontsize=12;

xxlim=[1e-3 1e3];
xxtick=power(10,-3:1:3);
xp=1.5e-3;

ReviseLineStyle='-';
ReviseColor='c';

figure;
BottomRetainWidth=0.05;
LeftRetainWidth=0.1;
Height=0.23;
Width=0.39;

% BottomRetainWidth=0.05;
% LeftRetainWidth=0.07;
% Height=0.23;
% Width=0.25;

%% lim
if strcmp(PluName,'1999CE119')   
%     yylimAI=[1e-5 1e3];
    yylimAI=[1e-7 1e1];
else
    yylimAI=[1e-7 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimAe=[1e-7 1e1];
else
    yylimAe=[1e-8 1e1];
end

if strcmp(PluName,'1999CE119')   
%     yylimGI=[1e-5 1e3];
    yylimGI=[1e-7 1e1];
else
    yylimGI=[1e-7 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimGe=[1e-7 1e1];
else
    yylimGe=[1e-8 1e1];
end

if strcmp(PluName,'1999CE119')   
%     yylimMI=[1e-10 1e5];
    yylimMI=[1e-13 1e1];
else
    yylimMI=[1e-12 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimMe=[1e-13 1e1];
else
    yylimMe=[1e-15 1e1];
end

labelx='$\nu_P$';


ReviseNo=length(RatioList1)+1;

% set(gcf,'Position',[400,100,700/4/0.618*3,700],'color','w');
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');

row=4;
%col=3;
col=2;

col1=[col+1 2*col+1 3*col+1]';
col2=col1+1;
col3=col2+1;

subplot(row,col,1);
TimesReviseList=TimesList(ReviseNo:end,1)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,TimesList,'k.');hold all;
loglog(RatioList2,TimesReviseList,'c.');
for ip=1:length(RatioList2)
    plot([RatioList2(ip) RatioList2(ip)],[TimesReviseList(ip) TimesList(ReviseNo+ip-1,1)],'c-');
end
loglog(RatioList1,Timesfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Timesfit2,'r-','linewidth',linewidth);
Timesfit1_extend=exp(polyval(PT1,log(RatioList2)));
loglog(RatioList2,Timesfit1_extend,'r-.','linewidth',linewidth);
if strcmp(PluName,'1999CE119')   
    yylim=[10 1e6];
else
    yylim=[1 1e5];
end
ylim(yylim);
xlim(xxlim);
plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

set(text(xp,5/10*yylim(2),['$$\ln{y} = ',num2str(PT1(1),'%.4f'),'\,\ln{x}+',num2str(PT1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,2/10*yylim(2),['$$R^2 = ',num2str(RT1(1,2)*RT1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(N)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% Res time
subplot(row,col,2);
loglog(RatioList,EjectTime,'k.');hold all;
% loglog(RatioList,MaxTime(:,1),'r.-');
% loglog(RatioList,MaxTime(:,2),'b.-');
% loglog(RatioList,MaxTime(:,3),'g.-');
yylim=[1e5 1e10];
ylim(yylim);
xlim(xxlim);
xxlim=get(gca,'xlim');
plot([barrier barrier],yylim,'r--');
ylabel('$T_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'xTick',[]);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% Ai
subplot(row,col,col1(1));
% SumReviseList=SumList(ReviseNo:end,1)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,SumList(:,1),'k.');hold all;
%loglog(RatioList,absDelList(:,1),'b.')

% loglog(RatioList2,SumReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[SumReviseList(ip) SumList(ReviseNo+ip-1,1)],'c-');
% end
loglog(RatioList1,Sumfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Sumfit2,'r-','linewidth',linewidth);
Sumfit1_extend=exp(polyval(PSum1,log(RatioList2)));
loglog(RatioList2,Sumfit1_extend,'r-.','linewidth',linewidth);

% ylabel('$A_I=\sum \Delta I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$A_I$','fontsize',fontsize,'Interpreter','latex');

yylim=yylimAI;
ylim(yylim);
xlim(xxlim);

plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

% set(gca,'yTick',power(10,-6:1:2));
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSum1(1),'%.4f'),'\,\ln{x}',num2str(PSum1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RSum1(1,2)*RSum1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
%% cumsum i
subplot(row,col,col1(2));
% CumsumReviseList=MaxCumsumList(ReviseNo:end,1)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,MaxCumsumList(:,1),'k.');hold all;
% loglog(RatioList,absDelList(:,1),'b.')
% loglog(RatioList,MaxList(:,1),'g.');

% loglog(RatioList2,CumsumReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[CumsumReviseList(ip) MaxCumsumList(ReviseNo+ip-1,1)],'c-');
% end
loglog(RatioList1,Cumsumfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Cumsumfit2,'r-','linewidth',linewidth);
Cumsumfit1_extend=exp(polyval(PCumsum1,log(RatioList2)));
loglog(RatioList2,Cumsumfit1_extend,'r-.','linewidth',linewidth);

% ylabel('$G_I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$G_I$','fontsize',fontsize,'Interpreter','latex');


yylim=yylimGI;
ylim(yylim);
xlim(xxlim);

plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

% set(gca,'yTick',power(10,-6:1:2));
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsum1(1),'%.4f'),'\,\ln{x}',num2str(PCumsum1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsum1(1,2)*RCumsum1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(B1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% M2
subplot(row,col,col1(3));
% RReviseList=RList(ReviseNo:end,1)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,RList(:,1),'k.');hold all;
% loglog(RatioList2,RReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[RReviseList(ip) RList(ReviseNo+ip-1,1)],'c-');
% end
loglog(RatioList1,Rfit1,'r-','linewidth',linewidth);
loglog(RatioList2,Rfit2,'r-','linewidth',linewidth);
Rfit1_extend=exp(polyval(PR1,log(RatioList2)));
loglog(RatioList2,Rfit1_extend,'r-.','linewidth',linewidth);

yylim=yylimMI;
ylim(yylim);
xlim(xxlim);
plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1)):2:log10(yylim(2))-1));

set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PR1(1),'%.4f'),'\,\ln{x}',num2str(PR1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,3/1000*yylim(2),['$$R^2 = ',num2str(RR1(1,2)*RR1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');

xlabel(labelx,'fontsize',fontsize,'Interpreter','latex');
% ylabel('$M_I=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_I$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025 0.05 0.05],'edgecolor','none','string',...,
           '(C1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(row,col,col2(1));
% SumReviseList=SumList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,SumList(:,2),'k.');hold all;
%loglog(RatioList,absDelList(:,2),'b.')

% loglog(RatioList2,SumReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[SumReviseList(ip) SumList(ReviseNo+ip-1,2)],'c-');
% end
loglog(RatioList1,Sumfitde1,'r-','linewidth',linewidth);
loglog(RatioList2,Sumfitde2,'r-','linewidth',linewidth);
Sumfitde1_extend=exp(polyval(PSumde1,log(RatioList2)));
loglog(RatioList2,Sumfitde1_extend,'r-.','linewidth',linewidth);

% ylabel('$A_e=\sum \Delta e$','fontsize',fontsize,'Interpreter','latex');
ylabel('$A_e$','fontsize',fontsize,'Interpreter','latex');


yylim=yylimAe;

ylim(yylim);
xlim(xxlim);

plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));
set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSumde1(1),'%.4f'),'\,\ln{x}',num2str(PSumde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RSumde1(1,2)*RSumde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
%%%%%% Cumsum
subplot(row,col,col2(2));
% MaxCumsumReviseList=MaxCumsumList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,MaxCumsumList(:,2),'k.');hold all;
% loglog(RatioList,absDelList(:,2),'b.')
% loglog(RatioList,MaxList(:,2),'g.');

% loglog(RatioList2,MaxCumsumReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[MaxCumsumReviseList(ip) MaxCumsumList(ReviseNo+ip-1,2)],'c-');
% end
loglog(RatioList1,Cumsumfitde1,'r-','linewidth',linewidth);
loglog(RatioList2,Cumsumfitde2,'r-','linewidth',linewidth);
Cumsumfitde1_extend=exp(polyval(PCumsumde1,log(RatioList2)));
loglog(RatioList2,Cumsumfitde1_extend,'r-.','linewidth',linewidth);

ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');


yylim=yylimGe;

ylim(yylim);
xlim(xxlim);

plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));
set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsumde1(1),'%.4f'),'\,\ln{x}',num2str(PCumsumde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsumde1(1,2)*RCumsumde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(B2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

subplot(row,col,col2(3));
% RReviseList=RList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,RList(:,2),'k.');hold all;
% loglog(RatioList2,RReviseList,'c.');
% for ip=1:length(RatioList2)
%     plot([RatioList2(ip) RatioList2(ip)],[RReviseList(ip) RList(ReviseNo+ip-1,2)],'c-');
% end
loglog(RatioList1,Rfitde1,'r-','linewidth',linewidth);
loglog(RatioList2,Rfitde2,'r-','linewidth',linewidth);
Rfitde1_extend=exp(polyval(PRde1,log(RatioList2)));
loglog(RatioList2,Rfitde1_extend,'r-.','linewidth',linewidth);



yylim=yylimMe;

ylim(yylim);
xlim(xxlim);
plot([barrier barrier],yylim,'r--');
hold off;
set(gca,'yTick',power(10,log10(yylim(1)):2:log10(yylim(2))-1));

set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PRde1(1),'%.4f'),'\,\ln{x}',num2str(PRde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,3/1000*yylim(2),['$$R^2 = ',num2str(RRde1(1,2)*RRde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');

% ylabel('$M_e=\sum {(\Delta e-\overline{\Delta e})}^{2}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_e$','fontsize',fontsize,'Interpreter','latex');

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025 0.05 0.05],'edgecolor','none','string',...,
           '(C2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
set(gca,'xTick',xxtick);
xlabel(labelx,'fontsize',fontsize,'Interpreter','latex');

%% Cumsum



% %%%%%%%%%%%%%%%%%%%%%%%%%%%da%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(row,col,col3(1));
% % SumReviseList=SumList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
% loglog(RatioList,SumList(:,3),'k.');hold all;
% %loglog(RatioList,absDelList(:,2),'b.')
% 
% % loglog(RatioList2,SumReviseList,'c.');
% % for ip=1:length(RatioList2)
% %     plot([RatioList2(ip) RatioList2(ip)],[SumReviseList(ip) SumList(ReviseNo+ip-1,2)],'c-');
% % end
% loglog(RatioList1,Sumfitda1,'r-','linewidth',linewidth);
% loglog(RatioList2,Sumfitda2,'r-','linewidth',linewidth);
% Sumfitda1_extend=exp(polyval(PSumda1,log(RatioList2)));
% loglog(RatioList2,Sumfitda1_extend,'r-.','linewidth',linewidth);
% 
% ylabel('$A_e=\sum \Delta e$','fontsize',fontsize,'Interpreter','latex');
% if strcmp(PluName,'1999CE119')   
%     yylim=[1e-7 1e1];
% else
%     yylim=[1e-8 1e1];
% end
% ylim(yylim);
% xlim(xxlim);
% 
% plot([barrier barrier],yylim,'r--');
% hold off;
% set(gca,'xTick',xxtick);
% 
% set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));
% set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSumda1(1),'%.4f'),'\,\ln{x}',num2str(PSumda1(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RSumda1(1,2)*RSumda1(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% set(gca,'xticklabel',[]);
% set(gca,'position',[3*LeftRetainWidth+2*Width BottomRetainWidth+2*Height Width Height]);
% annotation('textbox',[3*Width+LeftRetainWidth*1.5 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
%            '(A2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
% 
%        
% %%%%%% Cumsum
% subplot(row,col,col3(2));
% % MaxCumsumReviseList=MaxCumsumList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
% loglog(RatioList,MaxCumsumList(:,3),'k.');hold all;
% loglog(RatioList,absDelList(:,3),'b.')
% loglog(RatioList,MaxList(:,3),'g.');
% 
% % loglog(RatioList2,MaxCumsumReviseList,'c.');
% % for ip=1:length(RatioList2)
% %     plot([RatioList2(ip) RatioList2(ip)],[MaxCumsumReviseList(ip) MaxCumsumList(ReviseNo+ip-1,2)],'c-');
% % end
% loglog(RatioList1,Cumsumfitda1,'r-','linewidth',linewidth);
% loglog(RatioList2,Cumsumfitda2,'r-','linewidth',linewidth);
% Cumsumfitda1_extend=exp(polyval(PCumsumda1,log(RatioList2)));
% loglog(RatioList2,Cumsumfitda1_extend,'r-.','linewidth',linewidth);
% 
% ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
% 
% if strcmp(PluName,'1999CE119')   
%     yylim=[1e-7 1e1];
% else
%     yylim=[1e-8 1e1];
% end
% ylim(yylim);
% xlim(xxlim);
% 
% plot([barrier barrier],yylim,'r--');
% hold off;
% set(gca,'xTick',xxtick);
% 
% set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));
% set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsumda1(1),'%.4f'),'\,\ln{x}',num2str(PCumsumda1(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsumda1(1,2)*RCumsumda1(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% set(gca,'xticklabel',[]);
% set(gca,'position',[3*LeftRetainWidth+2*Width BottomRetainWidth+Height Width Height]);
% annotation('textbox',[3*Width+LeftRetainWidth*1.5 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
%            '(B2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
% 
% subplot(row,col,col3(3));
% % RReviseList=RList(ReviseNo:end,2)./EjectTime(ReviseNo:end)*1e9;
% loglog(RatioList,RList(:,3),'k.');hold all;
% % loglog(RatioList2,RReviseList,'c.');
% % for ip=1:length(RatioList2)
% %     plot([RatioList2(ip) RatioList2(ip)],[RReviseList(ip) RList(ReviseNo+ip-1,2)],'c-');
% % end
% loglog(RatioList1,Rfitda1,'r-','linewidth',linewidth);
% loglog(RatioList2,Rfitda2,'r-','linewidth',linewidth);
% Rfitda1_extend=exp(polyval(PRda1,log(RatioList2)));
% loglog(RatioList2,Rfitda1_extend,'r-.','linewidth',linewidth);
% 
% if strcmp(PluName,'1999CE119')   
%     yylim=[1e-13 1e1];
% else
%     yylim=[1e-15 1e1];
% end
% ylim(yylim);
% xlim(xxlim);
% plot([barrier barrier],yylim,'r--');
% hold off;
% set(gca,'yTick',power(10,log10(yylim(1)):2:log10(yylim(2))-1));
% 
% set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PRda1(1),'%.4f'),'\,\ln{x}',num2str(PRda1(2),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% set(text(xp,3/1000*yylim(2),['$$R^2 = ',num2str(RRda1(1,2)*RRda1(2,1),'%.4f'),'$$']),...,
%     'Interpreter','latex','fontsize',fontsize,'color','red');
% 
% ylabel('$M_2=\sum {(\Delta e-\overline{\Delta e})}^{2}$','fontsize',fontsize,'Interpreter','latex');
% 
% set(gca,'position',[3*LeftRetainWidth+2*Width BottomRetainWidth Width Height]);
% annotation('textbox',[3*Width+LeftRetainWidth*1.5 0.025 0.05 0.05],'edgecolor','none','string',...,
%            '(C2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
% set(gca,'xTick',xxtick);
% xlabel('$m/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
% 
% %% Cumsum
% 
