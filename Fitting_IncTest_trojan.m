%%Inc test
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end
% tag='data';

if ~strcmp(tag,'repeat')

clear;

Dir='ServerMount'; 

Npl=2;
DDir1='IncTestTro1';
DDir2='IncTestTro2';
PluName1='2004UP10';
PluName2='2006RJ103';

di_name='di_record_inout';
de_name='de_record_inout';

IncList=0:1:30;
IncList=IncList';

for ipl=1:Npl
    
    TimesList=zeros(length(IncList),1);
    SumList=zeros(length(IncList),2);
    StdList=zeros(length(IncList),2);
    RList=zeros(length(IncList),2);
    EndTime=zeros(length(IncList),1);
    plInc=zeros(length(IncList),1);
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
for i=1:length(IncList)
    Inc=IncList(i);
    disp(Inc);
    fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
    EndTime(i)=CE_record(end,1)/365;

    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
    plInc(i)=mean(plel(:,4));
    TimesList(i)=length(di_record_inout);
    
    SumList(i,1)=abs(sum(di_record_inout));
    StdList(i,1)=std(di_record_inout);
    RList(i,1)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumList(i,2)=abs(sum(de_record_inout));
    StdList(i,2)=std(de_record_inout);
    RList(i,2)=sum((de_record_inout-mean(de_record_inout)).^2);
end

    %sinInc0=1./sind(mean(tpInc));
    eval(['plIncMeanData',num2str(ipl),'=mean(plInc);']);
    eval(['TimesList',num2str(ipl),'=TimesList;']);
    eval(['SumList',num2str(ipl),'=SumList;']);
    eval(['Std',num2str(ipl),'=StdList;']);
    eval(['RList',num2str(ipl),'=RList;']);
%     Last=1;
%     Fitx=1./abs(sind(IncList(Last:end)-mean(tpInc)));
%     TimesListFit=TimesList(Last:end);
%     SumListFit=SumList(Last:end);
%     StdListFit=StdList(Last:end);
%     RListFit=RList(Last:end);

end

%plotx=1./abs(sind(IncList-mean(tpInc)));

end

A=(3/30-1/40)/2;
B=((2/30-1/40)/30)^(1/2);

%func=@(x)1./cosd(x);
%func=@(x)1./(A-B*cosd(x)).^(1/2).*sind(90-2*x);
func=@(x)sind(x);
plotf='semilogy';

plotx=func(IncList);
plIncMean1=func(plIncMeanData1);
plIncMean2=func(plIncMeanData2);


%xtick=[30 20 10 5 3 2 1 0.5 0.3 0.2 0.1];

markersize=15;
fontsize=12;

figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
row=4;
col=2;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.40;

subplot(row,col,1);
eval([plotf,'(1,1,''w'');']);hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
legend([h1 h2],{'1999CE119&2004UP10','1999CE119&2006RJ103'},'fontsize',fontsize,'location','northeast');

for ipl=1:Npl
    
    TimesList=eval(['TimesList',num2str(ipl)]);
    plIncMean=eval(['plIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx,TimesList(:,1),[color,''.-''],''markersize'',markersize);']);
    yylim=[1e2 1e4];
    ylim(yylim);

    plot([plIncMean plIncMean],yylim,'k--','linewidth',2);
    
end
hold off;
%semilogy(plotx(2:end),1./sind(plotx(2:end)-mean(tpInc))*2e2,'r-');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
%loglog(Fitx,Timesfit,'r-');
%set(gca,'xtick',1./sind(xtick));
%set(gca,'xticklabel',xtick);
%plot([sinInc0 sinInc0],yylim,'k--');
xlabel('$\sin{I_T}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'ytick',power(10,3:1:4));
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,3);
eval([plotf,'(1,3,''w'');']);hold all;
for ipl=1:Npl
    
    SumList=eval(['SumList',num2str(ipl)]);
    plIncMean=eval(['plIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx,SumList(:,1),[color,''.-''],''markersize'',markersize);']);
    yylim=[1e-6 1e0];
    ylim(yylim);
    plot([plIncMean plIncMean],yylim,'k--','linewidth',2);
    
%     [PSum,HSum]=polyfit(log(plotx),log(SumList(:,1)),1);  
%     RSum=corrcoef(log(plotx),log(SumList(:,1)));
%     Sumfit=exp(polyval(PSum,log(plotx)));
    
%     [PSum,HSum]=polyfit(plotx,log(SumList(:,1)),1);  
%     RSum=corrcoef(plotx,log(SumList(:,1)));
%     Sumfit=exp(polyval(PSum,plotx));
%     
%     eval([plotf,'(plotx,Sumfit,[color,''-''],''markersize'',markersize);']);

end
hold off;
%plot(plotx,SumList,'k.-','markersize',markersize);hold on;
%yylim=get(gca,'ylim');
%loglog(IncList,Sumfit,'r.-');
% 
% subplot(row,1,4);
% loglog(IncList,StdList,'k.-');hold on;
% loglog(IncList,Stdfit,'r.-');
xlabel('$\sin{I_T}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_I|=\left|\sum \Delta I\right|~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xticklabel',[]);
set(gca,'ytick',power(10,-5:1:-1));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,5);
eval([plotf,'(1,1,''w'');']);hold all;
for ipl=1:Npl
    
    RList=eval(['RList',num2str(ipl)]);
    plIncMean=eval(['plIncMean',num2str(ipl)]);

    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx,RList(:,1),[color,''.-''],''markersize'',markersize);']);
    yylim=[1e-7 1e-1];
    ylim(yylim);
    plot([plIncMean plIncMean],yylim,'k--','linewidth',2);
    
%     [PR,HR]=polyfit(log(plotx),log(RList(:,1)),1);  
%     RR=corrcoef(log(plotx),log(RList(:,1)));
%     Rfit=exp(polyval(PR,log(plotx)));
    
    [PR,HR]=polyfit(plotx,log(RList(:,1)),1);  
    RR=corrcoef(plotx,log(RList(:,1)));
    Rfit=exp(polyval(PR,plotx));
    
    disp(['ipl:',num2str(ipl)]);
    disp(['Corr:',num2str(RR(1,2)*RR(2,1))]);
    disp(['P1:',num2str(PR(1)),'P2',num2str(PR(2))]);
    
    eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
    
end
hold off;
xlabel('$\sin{I_T}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
%loglog(Fitx,Rfit,'r.-');
% set(gca,'xtick',1./sind(xtick));
% set(gca,'xticklabel',xtick);
%yylim=get(gca,'ylim');
%plot([sinInc0 sinInc0],yylim,'k--');
set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
set(gca,'ytick',power(10,-7:1:-2));
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+1*Height 0.03 0.03],'edgecolor','none','string',...,
           '(C1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
       
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
subplot(row,col,4);
eval([plotf,'(1,3,''w'');']);hold all;
for ipl=1:Npl
    
    SumList=eval(['SumList',num2str(ipl)]);
    plIncMean=eval(['plIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx,SumList(:,2),[color,''.-''],''markersize'',markersize);']);
    yylim=[1e-6 1e0];
    ylim(yylim);
    plot([plIncMean plIncMean],yylim,'k--','linewidth',2);
    
%     [PSume,HSume]=polyfit(plotx,log(SumList(:,2)),1);  
%     RSume=corrcoef(plotx,log(SumList(:,2)));
%     Sumfit=exp(polyval(PSume,plotx));
% 
%     eval([plotf,'(plotx,Sumfit,[color,''-''],''markersize'',markersize);']);

end
hold off;
%plot(plotx,SumList,'k.-','markersize',markersize);hold on;
%yylim=get(gca,'ylim');
%loglog(IncList,Sumfit,'r.-');
% 
% subplot(row,1,4);
% loglog(IncList,StdList,'k.-');hold on;
% loglog(IncList,Stdfit,'r.-');
xlabel('$\sin{I_T}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_e|=\left|\sum \Delta e\right|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xticklabel',[]);
set(gca,'ytick',power(10,-5:1:0));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,6);
eval([plotf,'(1,1,''w'');']);hold all;
for ipl=1:Npl
    
    RList=eval(['RList',num2str(ipl)]);
    plIncMean=eval(['plIncMean',num2str(ipl)]);

    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    % loglog(plotx,RList,'k.-');hold on;
    eval([plotf,'(plotx,RList(:,2),[color,''.-''],''markersize'',markersize);']);
    yylim=[1e-8 1e-4];
    ylim(yylim);
    
    plot([plIncMean plIncMean],yylim,'k--','linewidth',2);

    [PRe,HRe]=polyfit(plotx,log(RList(:,2)),1);  
    RRe=corrcoef(plotx,log(RList(:,2)));
    Rfit=exp(polyval(PRe,plotx));
    
    eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
    
    disp(['ipl:',num2str(ipl)]);
    disp(['Corr:',num2str(RRe(1,2)*RRe(2,1))]);
    disp(['P1:',num2str(PRe(1)),'P2',num2str(PRe(2))]);
end
hold off;
set(gca,'ytick',power(10,-8:1:-5));
xlabel('$\sin{I_T}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta e-\overline{\Delta e})}^{2}$','fontsize',fontsize,'Interpreter','latex');
%loglog(Fitx,Rfit,'r.-');
% set(gca,'xtick',1./sind(xtick));
% set(gca,'xticklabel',xtick);
%yylim=get(gca,'ylim');
%plot([sinInc0 sinInc0],yylim,'k--');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+Height 0.03 0.03],'edgecolor','none','string',...,
           '(C2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

