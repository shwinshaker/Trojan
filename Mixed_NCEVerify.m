%% IncPlutino
if ~exist('TimesList','var')

clear;
Dir='ServerMount'; 

Npl=2;

DDir1='IncTest';
DDir2='IncTest2';
PluName1='1999CE119';
PluName2='2001FU172';
di_name='di_record_inout';

IncList=0:1:30;
IncList=IncList';

aPL=zeros(length(IncList),2);
aTL=zeros(length(IncList),2);
ePL=zeros(length(IncList),2);
eTL=zeros(length(IncList),2);
IPL=zeros(length(IncList),2);
ITL=zeros(length(IncList),2);
TimesList=zeros(length(IncList),2);
tpIncMean=zeros(1,2);

for ipl=1:Npl
    
    tpInc=zeros(length(IncList),1);
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
    for i=1:length(IncList)
        Inc=IncList(i);
        disp(Inc);
        fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
        di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
        % CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
        
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
        tpInc(i)=mean(tpel(:,4));
        
        plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
        aPL(i,ipl)=mean(plel(:,2));
        aTL(i,ipl)=mean(tpel(:,2));
        ePL(i,ipl)=mean(plel(:,3));
        eTL(i,ipl)=mean(tpel(:,3));
        IPL(i,ipl)=mean(plel(:,4))/180*pi;
        ITL(i,ipl)=mean(tpel(:,4))/180*pi;
        
        TimesList(i,ipl)=length(di_record_inout);
    end
    
    %% sort IP
    [IPL(:,ipl),ind]=sort(IPL(:,ipl));
    TimesList(:,ipl)=TimesList(ind,ipl);
    aPL(:,ipl)=aPL(ind,ipl);
    aTL(:,ipl)=aTL(ind,ipl);
    ePL(:,ipl)=ePL(ind,ipl);
    eTL(:,ipl)=eTL(ind,ipl);
    ITL(:,ipl)=ITL(ind,ipl);
    
    tpIncMean(1,ipl)=mean(tpInc);
end
end

%% IncTrojan
if ~exist('TimesList_IT','var')

Dir='ServerMount'; 

Npl=2;

DDir1='IncTestTro1';
DDir2='IncTestTro2';
PluName1='2004UP10';
PluName2='2006RJ103';
di_name='di_record_inout';

IncList=0:1:30;
IncList=IncList';

aPL_IT=zeros(length(IncList),2);
aTL_IT=zeros(length(IncList),2);
ePL_IT=zeros(length(IncList),2);
eTL_IT=zeros(length(IncList),2);
IPL_IT=zeros(length(IncList),2);
ITL_IT=zeros(length(IncList),2);
TimesList_IT=zeros(length(IncList),2);
plIncMean_IT=zeros(1,2);

for ipl=1:Npl
    
    plInc=zeros(length(IncList),1);
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
    for i=1:length(IncList)
        Inc=IncList(i);
        disp(Inc);
        fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
        di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
        % CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
        
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
        plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
        plInc(i)=mean(plel(:,4));

        aPL_IT(i,ipl)=mean(plel(:,2));
        aTL_IT(i,ipl)=mean(tpel(:,2));
        ePL_IT(i,ipl)=mean(plel(:,3));
        eTL_IT(i,ipl)=mean(tpel(:,3));
        IPL_IT(i,ipl)=mean(plel(:,4))/180*pi;
        ITL_IT(i,ipl)=mean(tpel(:,4))/180*pi;
        
        TimesList_IT(i,ipl)=length(di_record_inout);
    end
    
    %% sort IP
    [ITL_IT(:,ipl),ind]=sort(ITL_IT(:,ipl));
    TimesList_IT(:,ipl)=TimesList_IT(ind,ipl);
    aPL_IT(:,ipl)=aPL_IT(ind,ipl);
    aTL_IT(:,ipl)=aTL_IT(ind,ipl);
    ePL_IT(:,ipl)=ePL_IT(ind,ipl);
    eTL_IT(:,ipl)=eTL_IT(ind,ipl);
    IPL_IT(:,ipl)=IPL_IT(ind,ipl);
    
    plIncMean_IT(1,ipl)=mean(plInc);
end
end

%% EccPlutino
if ~exist('TimesList_eP','var')

Dir='swiftdata/Trojan'; 

Npl=2;

DDir1='EccTest2';
DDir2='EccTest3';
PluName1='1999CE119';
PluName2='2001FU172';
di_name='di_record_inout';

EccList=0.2:0.005:0.4;
EccList=EccList';

aPL_eP=zeros(length(EccList),2);
aTL_eP=zeros(length(EccList),2);
ePL_eP=zeros(length(EccList),2);
eTL_eP=zeros(length(EccList),2);
IPL_eP=zeros(length(EccList),2);
ITL_eP=zeros(length(EccList),2);
TimesList_eP=zeros(length(EccList),2);

for ipl=1:Npl
    
    plInc=zeros(length(IncList),1);
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
    for i=1:length(IncList)
        Ecc=EccList(i);
        disp(Ecc);
        fname=[PluName,'_',num2str(sprintf('%.3f',Ecc)),'Ecc'];
        di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
        % CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
        
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
        plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
%         plInc(i)=mean(plel(:,4));

        aPL_eP(i,ipl)=mean(plel(:,2));
        aTL_eP(i,ipl)=mean(tpel(:,2));
        ePL_eP(i,ipl)=mean(plel(:,3));
        eTL_eP(i,ipl)=mean(tpel(:,3));
        IPL_eP(i,ipl)=mean(plel(:,4))/180*pi;
        ITL_eP(i,ipl)=mean(tpel(:,4))/180*pi;
        
        TimesList_eP(i,ipl)=length(di_record_inout);
    end
    
    %% sort IP
    [ePL_eP(:,ipl),ind]=sort(ePL_eP(:,ipl));
    TimesList_eP(:,ipl)=TimesList_eP(ind,ipl);
    aPL_eP(:,ipl)=aPL_eP(ind,ipl);
    aTL_eP(:,ipl)=aTL_eP(ind,ipl);
    eTL_eP(:,ipl)=eTL_eP(ind,ipl);
    IPL_eP(:,ipl)=IPL_eP(ind,ipl);
    ITL_eP(:,ipl)=ITL_eP(ind,ipl);
    
%     plIncMean_IT(1,ipl)=mean(plInc);
end
end

figure;

markersize=15;
fontsize=12;

Npl=2;

subplot(2,2,1);
semilogy(1,1,'w');hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
legend([h1 h2],{'1999CE119&2004UP10','2001FU172&2004UP10'},'fontsize',fontsize,'location','northeast');

NCE=zeros(length(IPL),2);
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    semilogy(IPL(:,ipl)/pi*180,TimesList(:,ipl),[color,'.-'],'markersize',markersize);
    yylim=[1e2 1e4];
    ylim(yylim);

    plot([tpIncMean(1,ipl) tpIncMean(1,ipl)],yylim,'k--','linewidth',2);
    NCE(:,ipl)=Fun_NCE(aPL(:,ipl),aTL(:,ipl),ePL(:,ipl),eTL(:,ipl),IPL(:,ipl),ITL(:,ipl),1);
    semilogy(IPL(:,ipl)/pi*180,NCE(:,ipl),[color,'--']);
    
end
hold off;
% semilogy(plotx(2:end),1./sind(plotx(2:end)-mean(tpInc))*2e2,'r-');
% set(gca,'xticklabel',[]);
% set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
%loglog(Fitx,Timesfit,'r-');
%set(gca,'xtick',1./sind(xtick));
%set(gca,'xticklabel',xtick);
%plot([sinInc0 sinInc0],yylim,'k--');
xlabel('$I_P~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'ytick',power(10,3:1:4));
%annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+3*Height 0.03 0.03],'edgecolor','none','string',...,
%            '(A)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(2,2,2);
semilogy(1,1,'w');hold all;
NCE_IT=zeros(length(ITL_IT),2);
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    semilogy(ITL_IT(:,ipl)/pi*180,TimesList_IT(:,ipl),[color,'.-'],'markersize',markersize);
    yylim=[1e2 1e4];
    ylim(yylim);

    plot([plIncMean_IT(1,ipl) plIncMean_IT(1,ipl)],yylim,'k--','linewidth',2);
    NCE_IT(:,ipl)=Fun_NCE(aPL_IT(:,ipl),aTL_IT(:,ipl),ePL_IT(:,ipl),eTL_IT(:,ipl),IPL_IT(:,ipl),ITL_IT(:,ipl),1);
    semilogy(ITL_IT(:,ipl)/pi*180,NCE_IT(:,ipl),[color,'--']);
    
end
hold off;
xlabel('$I_T~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

subplot(2,2,3);
semilogy(1,1,'w');hold all;
NCE_eP=zeros(length(ePL_eP),2);
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    semilogy(ePL_eP(:,ipl),TimesList_eP(:,ipl),[color,'.-'],'markersize',markersize);
    yylim=[1e2 1e4];
    ylim(yylim);

%     plot([plIncMean_IT(1,ipl) plIncMean_IT(1,ipl)],yylim,'k--','linewidth',2);
    NCE_eP(:,ipl)=Fun_NCE(aPL_eP(:,ipl),aTL_eP(:,ipl),ePL_eP(:,ipl),eTL_eP(:,ipl),IPL_eP(:,ipl),ITL_eP(:,ipl),1);
    semilogy(ePL_eP(:,ipl),NCE_eP(:,ipl),[color,'--']);
    
end
hold off;
xlabel('$I_T~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

