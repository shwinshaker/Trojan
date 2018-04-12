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
DDir1='IncTest';
DDir2='IncTest2';
PluName1='1999CE119';
PluName2='2001FU172';

di_name='di_record_inout';
de_name='de_record_inout';
%IncList=1./exp(log(1):0.1:log(1000));
IncList=0:1:30;
IncList=IncList';

aPL=zeros(length(IncList),2);
aTL=zeros(length(IncList),2);
eTL=zeros(length(IncList),2);
ePL=zeros(length(IncList),2);
IPL=zeros(length(IncList),2);
ITL=zeros(length(IncList),2);

TimesList=zeros(length(IncList),2);
SumIList=zeros(length(IncList),2);
StdIList=zeros(length(IncList),2);
RIList=zeros(length(IncList),2);

SumeList=zeros(length(IncList),2);
StdeList=zeros(length(IncList),2);
ReList=zeros(length(IncList),2);

EjectTime=zeros(length(IncList),2);
EndTime=zeros(length(IncList),2);
tpInc=zeros(length(IncList),2);
    
for ipl=1:Npl
   
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
for i=1:length(IncList)
    Inc=IncList(i);
    disp(Inc);
    fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
    
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
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

    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');

    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);
    di_record_inout=di_record_inout(1:ejectCENo)/180*pi;
    de_record_inout=de_record_inout(1:ejectCENo);
    
    CEtime=CE_record(1:ejectCENo,1)/365.25;
    EndTime(i,ipl)=CEtime(end);
    
    tpInc(i,ipl)=mean(tpel(:,4));
    
    aPL(i,ipl)=mean(plel(:,2));
    aTL(i,ipl)=mean(tpel(:,2));
    ePL(i,ipl)=mean(plel(:,3));
    eTL(i,ipl)=mean(tpel(:,3));
    IPL(i,ipl)=mean(plel(:,4))/180*pi;
    ITL(i,ipl)=mean(tpel(:,4))/180*pi;
    
    EjectTime(i,ipl)=ejecttime/365.25;
    TimesList(i,ipl)=length(di_record_inout);
    
    SumIList(i,ipl)=abs(sum(di_record_inout));
    StdIList(i,ipl)=std(di_record_inout);
    RIList(i,ipl)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumeList(i,ipl)=abs(sum(de_record_inout));
    StdeList(i,ipl)=std(de_record_inout);
    ReList(i,ipl)=sum((de_record_inout-mean(de_record_inout)).^2);
end

    %sinInc0=1./sind(mean(tpInc));
%     tpIncMean(ipl)=mean(tpInc(:,ipl));
%     eval(['TimesList',num2str(ipl),'=TimesList;']);
%     eval(['SumList',num2str(ipl),'=SumList;']);
%     eval(['Std',num2str(ipl),'=StdList;']);
%     eval(['RList',num2str(ipl),'=RList;']);
   
    %% sort IP
    [IPL(:,ipl),ind]=sort(IPL(:,ipl));
    
    aPL(:,ipl)=aPL(ind,ipl);
    aTL(:,ipl)=aTL(ind,ipl);
    ePL(:,ipl)=ePL(ind,ipl);
    eTL(:,ipl)=eTL(ind,ipl);
    ITL(:,ipl)=ITL(ind,ipl);
    
    TimesList(:,ipl)=TimesList(ind,ipl);
    EjectTime(:,ipl)=EjectTime(ind,ipl);
    SumIList(:,ipl)=SumIList(ind,ipl);
    RIList(:,ipl)=RIList(ind,ipl);
    StdIList(:,ipl)=StdIList(ind,ipl);
    SumeList(:,ipl)=SumeList(ind,ipl);
    ReList(:,ipl)=ReList(ind,ipl);
    StdeList(:,ipl)=StdeList(ind,ipl);

%     Last=1;
%     Fitx=1./abs(sind(IncList(Last:end)-mean(tpInc)));
%     TimesListFit=TimesList(Last:end);
%     SumListFit=SumList(Last:end);
%     StdListFit=StdList(Last:end);
%     RListFit=RList(Last:end);

end

end


%% Theo

if ~(exist('NCE','var') && exist('MI','var') && exist('Me','var'))
    
    %% Theo
    NCE=zeros(length(IPL),2);
    for ipl=1:Npl
        %NCE(:,ipl)=Fun_NCE(aPL(:,ipl),aTL(:,ipl),ePL(:,ipl),eTL(:,ipl),IPL(:,ipl),ITL(:,ipl),1);
        NCE(:,ipl)=Fun_CollisionProb(aTL(:,ipl),aPL(:,ipl),ePL(:,ipl),IPL(:,ipl))*1e9/10;
    end
    
    MI=zeros(length(NCE(:,ipl)),2);
    AI=zeros(length(NCE(:,ipl)),2);
    for ipl=1:Npl
        for id=1:length(NCE(:,ipl))
            [dinc,~]=Fun_diDstb_theo(100000,1,...,
                aPL(id,ipl),ePL(id,ipl),IPL(id,ipl)/pi*180,...,
                aTL(id,ipl),eTL(id,ipl),ITL(id,ipl)/pi*180);
            dinc=dinc/180*pi;
            MI(id,ipl)=std(dinc)^2*NCE(id,ipl);
        end
        AI(:,ipl)=sqrt(2/pi*MI(:,ipl));
    end

    Me=zeros(length(NCE(:,ipl)),2);
    Ae=zeros(length(NCE(:,ipl)),2);
    for ipl=1:Npl
        for id=1:length(NCE(:,ipl))
            [~,de]=Fun_diDstb_theo(100000,1,...,
                aPL(id,ipl),ePL(id,ipl),IPL(id,ipl)/pi*180,...,
                aTL(id,ipl),eTL(id,ipl),ITL(id,ipl)/pi*180);
            Me(id,ipl)=std(de)^2*NCE(id,ipl);
        end
        Ae(:,ipl)=sqrt(2/pi*Me(:,ipl));
    end
    
end

% A=(3/30-1/40)/2;
% B=((2/30-1/40)/30)^(1/2);

%func=@(x)1./cosd(x);
%func=@(x)1./(A-B*cosd(x)).^(1/2).*sind(90-2*x);
func=@(x)sind(x);
plotx=zeros(length(IPL(:,1)),2);
for ipl=1:Npl
    plotx(:,ipl)=func(IPL(:,ipl)/pi*180);
end
%func=@(x)x;
tpIncMean=func(mean(tpInc,1));

plotf='semilogy';

% % plotx=func(IncList);
% tpIncMean1=func(tpIncMeanData1);
% tpIncMean2=func(tpIncMeanData2);

%xtick=[30 20 10 5 3 2 1 0.5 0.3 0.2 0.1];

markersize=12;
linewidth=2;
fontsize=12;

figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
row=4;
col=2;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.40;

yylimA=[1e-8 1e-1];

subplot(row,col,1);
eval([plotf,'(0,1,''w'');']);hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
legend([h1 h2],{'1999CE119&2004UP10','2001FU172&2004UP10'},'fontsize',fontsize,'location','northeast');

for ipl=1:Npl
        
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx(:,ipl),TimesList(:,ipl),[color,''.''],''markersize'',markersize);']);
    
    semilogy(plotx(:,ipl),NCE(:,ipl),[color,'.-'],'linewidth',linewidth);
    
    yylim=[1e2 1e4];
    ylim(yylim);
    plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylim,'k--','linewidth',2);

end
hold off;
%semilogy(plotx(2:end),1./sind(plotx(2:end)-mean(tpInc))*2e2,'r-');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
%loglog(Fitx,Timesfit,'r-');
%set(gca,'xtick',1./sind(xtick));
%set(gca,'xticklabel',xtick);
%plot([sinInc0 sinInc0],yylim,'k--');
% xlabel('$\sin{I_P}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'ytick',power(10,3:1:4));
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(N)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

%% Res time
subplot(row,col,2);
eval([plotf,'(0,1,''w'');']);hold all;
for ipl=1:Npl
        
    switch ipl
        case 1
            color='r';
        case 2
            color='b';
    end
    semilogy(plotx(:,ipl),EjectTime(:,ipl),[color,'.'],'markersize',markersize);hold all;
    yylim=[1e6 1e10];
    ylim(yylim);
    plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylim,'k--','linewidth',2);

end


set(gca,'xticklabel',[]);
ylabel('$t_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
% set(gca,'xTick',[]);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

% subplot(row,col,5);
% eval([plotf,'(0,1e-5,''w'');']);hold all;
% for ipl=1:Npl
%     
%     switch ipl
%         case 1 
%             color='r';
%         case 2
%             color='b';
%     end
%     eval([plotf,'(plotx(:,ipl),RIList(:,ipl),[color,''.-''],''markersize'',markersize);']);
%     
%     yylim=get(gca,'ylim');
%     %yylim=[1e-7 1e-1];
%     ylim(yylim);
%     plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylim,'k--','linewidth',2);
%     
%     %% Theo
% %     aRL=aTL(:,ipl)./aPL(:,ipl);
% %     A=(3-aRL)/2;
% %     B=(2-aRL).^(1/2);
% %     kappa=3.5;
% %     mP=6.56e-9;
% %     ChiM=2*3^(-5/3)*(aRL/kappa).^2*mP^(4/3)./(A-B).*NCE(:,ipl)*1;
% %     MI=ChiM*1.*(-1+log(NCE(:,ipl)));
%     
%     semilogy(plotx(:,ipl),MI(:,ipl),[color,'--']);
%     
% %     %% Fit
% % %     [PR,HR]=polyfit(log(plotx),log(RList(:,1)),1);  
% % %     RR=corrcoef(log(plotx),log(RList(:,1)));
% % %     Rfit=exp(polyval(PR,log(plotx)));
% %     
% %     [PR,HR]=polyfit(plotx,log(RList(:,1)),1);  
% %     RR=corrcoef(plotx,log(RList(:,1)));
% %     Rfit=exp(polyval(PR,plotx));
% %     
% %     disp(['ipl:',num2str(ipl)]);
% %     disp(['Corr:',num2str(RR(1,2)*RR(2,1))]);
% %     disp(['P1:',num2str(PR(1)),'P2',num2str(PR(2))]);
% %     
% %     eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
%     
% end
% hold off;
% xlabel('$\sin{I_P}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$M_I$','fontsize',fontsize,'Interpreter','latex');
% %loglog(Fitx,Rfit,'r.-');
% % set(gca,'xtick',1./sind(xtick));
% % set(gca,'xticklabel',xtick);
% %yylim=get(gca,'ylim');
% %plot([sinInc0 sinInc0],yylim,'k--');
% set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
% %set(gca,'ytick',power(10,-7:1:-2));
% annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+1*Height 0.03 0.03],'edgecolor','none','string',...,
%            '(C1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
%        
       
subplot(row,col,3);
eval([plotf,'(0,1e-5,''w'');']);hold all;
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
%     eval([plotf,'(plotx(:,ipl),SumIList(:,ipl),[color,''.-''],''markersize'',markersize);']);
    eval([plotf,'(plotx(:,ipl),SumIList(:,ipl),[color,''.''],''markersize'',markersize);']);

    %     yylim=get(gca,'ylim');
    ylim(yylimA);
    plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylimA,'k--','linewidth',2);
    
    %% Theo %% Random Walk
    semilogy(plotx(:,ipl),AI(:,ipl),[color,'.-'],'linewidth',linewidth);

    
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
% xlabel('$\sin{I_P}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
xlabel('$\sin{I_P}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_I|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
% set(gca,'xticklabel',[]);
set(gca,'yTick',power(10,log10(yylimA(1))+1:log10(yylimA(2))-1));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');


       
%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
subplot(row,col,4);
eval([plotf,'(1,3,''w'');']);hold all;
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx(:,ipl),SumeList(:,ipl),[color,''.''],''markersize'',markersize);']);
%     eval([plotf,'(plotx(:,ipl),SumeList(:,ipl),[color,''.'']);']);

    ylim(yylimA);
    plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylimA,'k--','linewidth',2);
    
    semilogy(plotx(:,ipl),Ae(:,ipl),[color,'-'],'linewidth',linewidth);
       
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
xlabel('$\sin{I_P}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_e|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
% set(gca,'xticklabel',[]);
% set(gca,'ytick',power(10,-5:1:0));
set(gca,'yTick',power(10,log10(yylimA(1))+1:log10(yylimA(2))-1));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

% subplot(row,col,6);
% eval([plotf,'(1,1,''w'');']);hold all;
% for ipl=1:Npl
%     
% 
%     switch ipl
%         case 1 
%             color='r';
%         case 2
%             color='b';
%     end
%     % loglog(plotx,RList,'k.-');hold on;
%     eval([plotf,'(plotx(:,ipl),ReList(:,ipl),[color,''.-''],''markersize'',markersize);']);
%     yylim=[1e-8 1e-4];
%     ylim(yylim);
%     plot([tpIncMean(:,ipl) tpIncMean(:,ipl)],yylim,'k--','linewidth',2);
% 
%     semilogy(plotx(:,ipl),Me(:,ipl),[color,'--'])
% 
%     %% Fit
% %     [PRe,HRe]=polyfit(plotx,log(RList(:,2)),1);  
% %     RRe=corrcoef(plotx,log(RList(:,2)));
% %     Rfit=exp(polyval(PRe,plotx));
% %     
% %     eval([plotf,'(plotx,Rfit,[color,''-''],''markersize'',markersize);']);
% %     
% %     disp(['ipl:',num2str(ipl)]);
% %     disp(['Corr:',num2str(RRe(1,2)*RRe(2,1))]);
% %     disp(['P1:',num2str(PRe(1)),'P2',num2str(PRe(2))]);
% end
% hold off;
% set(gca,'ytick',power(10,-8:1:-5));
% xlabel('$\sin{I_P}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$M_e$','fontsize',fontsize,'Interpreter','latex');
% %loglog(Fitx,Rfit,'r.-');
% % set(gca,'xtick',1./sind(xtick));
% % set(gca,'xticklabel',xtick);
% %yylim=get(gca,'ylim');
% %plot([sinInc0 sinInc0],yylim,'k--');
% set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
% annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+Height 0.03 0.03],'edgecolor','none','string',...,
%            '(C2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
% 
