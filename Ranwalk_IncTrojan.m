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
plInc=zeros(length(IncList),2);

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

    plInc(i,ipl)=mean(plel(:,4));
    
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

    %% sort IT
    [ITL(:,ipl),ind]=sort(ITL(:,ipl));
    
    aPL(:,ipl)=aPL(ind,ipl);
    aTL(:,ipl)=aTL(ind,ipl);
    ePL(:,ipl)=ePL(ind,ipl);
    eTL(:,ipl)=eTL(ind,ipl);
    IPL(:,ipl)=IPL(ind,ipl);
    
    TimesList(:,ipl)=TimesList(ind,ipl);
    EjectTime(:,ipl)=EjectTime(ind,ipl);
    SumIList(:,ipl)=SumIList(ind,ipl);
    RIList(:,ipl)=RIList(ind,ipl);
    StdIList(:,ipl)=StdIList(ind,ipl);
    SumeList(:,ipl)=SumeList(ind,ipl);
    ReList(:,ipl)=ReList(ind,ipl);
    StdeList(:,ipl)=StdeList(ind,ipl);

end

%plotx=1./abs(sind(IncList-mean(tpInc)));

end

%% Theo

if ~(exist('NCE','var') && exist('MI','var') && exist('Me','var'))
    
    %% Theo
    NCE=zeros(length(ITL),2);
    for ipl=1:Npl
        NCE(:,ipl)=Fun_NCE(aPL(:,ipl),aTL(:,ipl),ePL(:,ipl),eTL(:,ipl),IPL(:,ipl),ITL(:,ipl),1);
    end
    
    MI=zeros(length(NCE(:,ipl)),2);
    AI=zeros(length(NCE(:,ipl)),2);
    AISim=zeros(length(NCE(:,ipl)),2);
    
    Me=zeros(length(NCE(:,ipl)),2);
    Ae=zeros(length(NCE(:,ipl)),2);
    
    for ipl=1:Npl
        for id=1:length(NCE(:,ipl))
            [dinc,de]=Fun_diDstb_theo(100000,1,...,
                aPL(id,ipl),ePL(id,ipl),IPL(id,ipl)/pi*180,...,
                aTL(id,ipl),eTL(id,ipl),ITL(id,ipl)/pi*180);
            dinc=dinc/180*pi;
            MI(id,ipl)=std(dinc)^2*NCE(id,ipl);
            Me(id,ipl)=std(de)^2*NCE(id,ipl);
            AISim(id,ipl)=Fun_totEff_sim(1,aPL(id,ipl),ePL(id,ipl),IPL(id,ipl),...,
                                           aTL(id,ipl),eTL(id,ipl),ITL(id,ipl));
        end
        AI(:,ipl)=sqrt(2/pi*MI(:,ipl));
        Ae(:,ipl)=sqrt(2/pi*Me(:,ipl));
    end
    AeSim=AISim;
 
end

plotf='semilogy';

func=@(x)sind(x);
plotx=func(ITL/pi*180);

plIncMean=func(mean(plInc,1));

markersize=12;
linewidth=2;
fontsize=12;

yylimAI=[1e-7 1e-1];

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
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx(:,ipl),TimesList(:,ipl),[color,''.''],''markersize'',markersize);']);
    
    semilogy(plotx(:,ipl),NCE(:,ipl),[color,'.-'],'linewidth',linewidth);
    
    yylim=[1e2 2e4];
    ylim(yylim);
    plot([plIncMean(:,ipl) plIncMean(:,ipl)],yylim,'k--','linewidth',2);
    
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
    plot([plIncMean(:,ipl) plIncMean(:,ipl)],yylim,'k--','linewidth',2);

end

set(gca,'xticklabel',[]);
ylabel('$t_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
% set(gca,'xTick',[]);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
subplot(row,col,3);
eval([plotf,'(0,1e-5,''w'');']);hold all;
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(plotx(:,ipl),SumIList(:,ipl),[color,''.''],''markersize'',markersize);']);
    
    semilogy(plotx(:,ipl),AI(:,ipl),[color,'.-'],'linewidth',linewidth);
    semilogy(plotx(:,ipl),AISim(:,ipl),'m--','linewidth',linewidth);

    yylim=yylimAI;
    ylim(yylim);
    plot([plIncMean(:,ipl) plIncMean(:,ipl)],yylim,'k--','linewidth',2);
    

end
hold off;
xlabel('$\sin{I_T}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_I|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
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
    semilogy(plotx(:,ipl),Ae(:,ipl),[color,'.-'],'linewidth',linewidth);
    semilogy(plotx(:,ipl),AeSim(:,ipl),'m--','linewidth',linewidth);
       
    yylim=yylimAI;
    ylim(yylim);
    plot([plIncMean(:,ipl) plIncMean(:,ipl)],yylim,'k--','linewidth',2);
    
end
hold off;
xlabel('$\sin{I_T}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_e|$','fontsize',fontsize,'Interpreter','latex');
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
