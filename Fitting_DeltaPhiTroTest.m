%%Inc test
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')

clear;

Dir='ServerMount'; 

Npl=2;
DDir1='DeltaPhiTestTro1';
DDir2='DeltaPhiTestTro2';
TroName1='2004UP10';
TroName2='2006RJ103';
%PluName2='2001FU172';

di_name='di_record_inout';
de_name='de_record_inout';


List=1:100;

for ipl=1:Npl
    
    TimesList=zeros(length(List),1);
    SumList=zeros(length(List),2);
    StdList=zeros(length(List),2);
    RList=zeros(length(List),2);
    EndTime=zeros(length(List),1);
    DeltaPhiList=zeros(length(List),1);
    
    DDir=eval(['DDir',num2str(ipl)]);
    TroName=eval(['TroName',num2str(ipl)]);
    
    errid=[];

for i=1:length(List)
    
    disp(['#',num2str(i)]);
    fname=[TroName,'_',num2str(i)];
    
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
    EndTime(i)=CE_record(end,1)/365;
    try
        plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
    catch
        errid=[errid;i];
        continue;
    end
    %% remove Plutinos exile from resonance
    temp=find(plel(:,2)>40.0 | plel(:,2)<39.0,1,'first');
    if ~isempty(temp)
        errid=[errid;i];
        continue;
    end

    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if ~isempty(temp)
        errid=[errid;i];
        continue;
    end
    
    Nepel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/Nepel.txt']);
    %phi=mod(3*(plel(:,5)+plel(:,6)+plel(:,7))-2*(Nepel(:,5)+Nepel(:,6)+Nepel(:,7))-(plel(:,5)+plel(:,6)),360);
    phi=mod((tpel(:,5)+tpel(:,6)+tpel(:,7))-(Nepel(:,5)+Nepel(:,6)+Nepel(:,7)),360);
    DeltaPhiList(i)=(max(phi)-min(phi))/2;
    
    TimesList(i)=length(di_record_inout);
    
    SumList(i,1)=abs(sum(di_record_inout));
    StdList(i,1)=std(di_record_inout);
    RList(i,1)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumList(i,2)=abs(sum(de_record_inout));
    StdList(i,2)=std(de_record_inout);
    RList(i,2)=sum((de_record_inout-mean(de_record_inout)).^2);
    
end

    %% remove err
    DeltaPhiList(errid)=[];
    TimesList(errid)=[];
    SumList(errid,:)=[];
    RList(errid,:)=[];
    StdList(errid,:)=[];
    
    %% sort DeltaPhi
    [DeltaPhiList,ind]=sort(DeltaPhiList);
    
    TimesList=TimesList(ind,:);
    SumList=SumList(ind,:);
    RList=RList(ind,:);
    StdList=StdList(ind,:);
    
    %% Add nan to separate tadpole and horseshoe
    Indsep=find(DeltaPhiList(2:end)-DeltaPhiList(1:end-1)>30);
    DeltaPhiList=[DeltaPhiList(1:Indsep);nan;DeltaPhiList(Indsep+1:end)];
    TimesList=[TimesList(1:Indsep);nan;TimesList(Indsep+1:end)];
    SumList=[SumList(1:Indsep,:);[nan nan];SumList(Indsep+1:end,:)];
    RList=[RList(1:Indsep,:);[nan nan];RList(Indsep+1:end,:)];
    StdList=[StdList(1:Indsep,:);[nan nan];StdList(Indsep+1:end,:)];
    
    % eval(['tpIncMeanData',num2str(ipl),'=mean(tpInc);']);
    eval(['TimesList',num2str(ipl),'=TimesList;']);
    eval(['SumList',num2str(ipl),'=SumList;']);
    eval(['Std',num2str(ipl),'=StdList;']);
    eval(['RList',num2str(ipl),'=RList;']);
    eval(['DeltaPhiList',num2str(ipl),'=DeltaPhiList;']);

end

end

func=@(x)sind(x);
plotf='semilogy';

% plotx=func(IncList);
% tpIncMean1=func(tpIncMeanData1);
% tpIncMean2=func(tpIncMeanData2);

%xtick=[30 20 10 5 3 2 1 0.5 0.3 0.2 0.1];

markersize=15;
fontsize=12;

yylim1=[1e2 1e4];
yylim4=[1e-7 1];
yylim6=[1e-8 1e-3];


%%  Figure -------------------------------------------------------
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
grid on;
legend([h1 h2],{'1999CE119&2004UP10','1999CE119&2006RJ103'},...,
    'fontsize',fontsize,'location','southwest');

for ipl=1:Npl
    
    DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);
    TimesList=eval(['TimesList',num2str(ipl)]);
    % tpIncMean=eval(['tpIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList,TimesList(:,1),[color,''.-''],''markersize'',markersize);']);
%     yylim=[1e2 1e4];
%     ylim(yylim);

%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
    
end
hold off;
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);

%set(gca,'xtick',1./sind(xtick));
%set(gca,'xticklabel',xtick);

xlim([0 180]);
set(gca,'xtick',0:30:180);
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
if exist('yylim1','var')
    yylim=yylim1;
else
    yylim=get(gca,'ylim');
end
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,3);
eval([plotf,'(1,3,''w'');']);hold all;
grid on;
for ipl=1:Npl
    
    DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);    
    SumList=eval(['SumList',num2str(ipl)]);
%     tpIncMean=eval(['tpIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList,SumList(:,1),[color,''.-''],''markersize'',markersize);']);
%     yylim=[1e-6 1e0];
%     ylim(yylim);
%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
    
end
hold off;

xlim([0 180]);
ylabel('$|A_I|=\left|\sum \Delta I\right|~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xtick',0:30:180);
set(gca,'xticklabel',[]);
yylim=get(gca,'ylim');
% yytick=get(gca,'yticklabel');
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,5);
eval([plotf,'(1,1,''w'');']);hold all;
grid on;
for ipl=1:Npl
    
    DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);        
    RList=eval(['RList',num2str(ipl)]);
%     tpIncMean=eval(['tpIncMean',num2str(ipl)]);

    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList,RList(:,1),[color,''.-''],''markersize'',markersize);']);
%     yylim=[1e-7 1e-1];
%     ylim(yylim);
%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
        
%     [PR,HR]=polyfit(plotx,log(RList(:,1)),1);  
%     RR=corrcoef(plotx,log(RList(:,1)));
%     Rfit=exp(polyval(PR,plotx));
%     
%     disp(['ipl:',num2str(ipl)]);
%     disp(['Corr:',num2str(RR(1,2)*RR(2,1))]);
%     disp(['P1:',num2str(PR(1)),'P2',num2str(PR(2))]);
%     
%     eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
    
end
hold off;

xlim([0 180]);
set(gca,'xtick',0:30:180);
xlabel('$\Delta\phi~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
yylim=get(gca,'ylim');
set(gca,'yTick',power(10,log10(yylim(1)):1:log10(yylim(2))-1));
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+1*Height 0.03 0.03],'edgecolor','none','string',...,
           '(C1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
       
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
subplot(row,col,4);
eval([plotf,'(1,3,''w'');']);hold all;
grid on;
for ipl=1:Npl
    
   DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);        
   SumList=eval(['SumList',num2str(ipl)]);
%     tpIncMean=eval(['tpIncMean',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList,SumList(:,2),[color,''.-''],''markersize'',markersize);']);
%     yylim=[1e-6 1e0];
%     ylim(yylim);
%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
    
end
hold off;

xlim([0 180]);
ylabel('$|A_e|=\left|\sum \Delta e\right|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xtick',0:30:180);
set(gca,'xticklabel',[]);
if exist('yylim4','var')
    yylim=yylim4;
else
    yylim=get(gca,'ylim');
end
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,6);
eval([plotf,'(1,1,''w'');']);hold all;
grid on;
for ipl=1:Npl
    
    DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);        
    RList=eval(['RList',num2str(ipl)]);
%     tpIncMean=eval(['tpIncMean',num2str(ipl)]);

    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    
    eval([plotf,'(DeltaPhiList,RList(:,2),[color,''.-''],''markersize'',markersize);']);
%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
%     yylim=[1e-8 1e-4];
%     ylim(yylim);
%     
%     [PRe,HRe]=polyfit(plotx,log(RList(:,2)),1);  
%     RRe=corrcoef(plotx,log(RList(:,2)));
%     Rfit=exp(polyval(PRe,plotx));
%     
%     eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
%     
%     disp(['ipl:',num2str(ipl)]);
%     disp(['Corr:',num2str(RRe(1,2)*RRe(2,1))]);
%     disp(['P1:',num2str(PRe(1)),'P2',num2str(PRe(2))]);
end
hold off;

xlim([0 180]);
set(gca,'xtick',0:30:180);
set(gca,'ytick',power(10,-8:1:-5));
xlabel('$\Delta\phi~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta e-\overline{\Delta e})}^{2}$','fontsize',fontsize,'Interpreter','latex');
if exist('yylim6','var')
    yylim=yylim6;
else
    yylim=get(gca,'ylim');
end
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1)):1:log10(yylim(2))-1));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+Height 0.03 0.03],'edgecolor','none','string',...,
           '(C2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

