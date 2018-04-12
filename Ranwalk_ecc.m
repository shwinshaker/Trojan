%%Ecc test
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')

clear;

Dir='swiftdata'; 
Npl=2;

DDir1='EccTest2';
DDir2='EccTest3';

PluName1='1999CE119';
PluName2='2001FU172';
di_name='di_record_inout';
de_name='de_record_inout';

EccList=0.2:0.005:0.4;
EccList=EccList';

aPL=zeros(length(EccList),2);
aTL=zeros(length(EccList),2);
ePL=zeros(length(EccList),2);
eTL=zeros(length(EccList),2);
IPL=zeros(length(EccList),2);
ITL=zeros(length(EccList),2);

EjectTime=zeros(length(EccList),2);
TimesList=zeros(length(EccList),2);
SumIList=zeros(length(EccList),2);
StdIList=zeros(length(EccList),2);
RIList=zeros(length(EccList),2);

SumeList=zeros(length(EccList),2);
StdeList=zeros(length(EccList),2);
ReList=zeros(length(EccList),2);

for ipl=1:Npl
        
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
for i=1:length(EccList)
    Ecc=EccList(i);
    disp(Ecc);

    fname=[PluName,'_',num2str(sprintf('%.3f',Ecc)),'Ecc'];

    tpel=load(['~/Documents/',Dir,'/Trojan/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);   
    plel=load(['~/Documents/',Dir,'/Trojan/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
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

    CE_record=load(['~/Documents/',Dir,'/Trojan/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
    if isempty(CE_record)
        di_record_inout=0;
        de_record_inout=0;
        TimesList(i,ipl)=0;
    else
        ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
        
        di_record_inout=load(['~/Documents/',Dir,'/Trojan/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
        de_record_inout=load(['~/Documents/',Dir,'/Trojan/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);
        di_record_inout=di_record_inout(1:ejectCENo)/180*pi;
        de_record_inout=de_record_inout(1:ejectCENo);
        
        TimesList(i,ipl)=length(di_record_inout);
    end
    
    aPL(i,ipl)=mean(plel(:,2));
    aTL(i,ipl)=mean(tpel(:,2));
    ePL(i,ipl)=mean(plel(:,3));
    eTL(i,ipl)=mean(tpel(:,3));
    IPL(i,ipl)=mean(plel(:,4))/180*pi;
    ITL(i,ipl)=mean(tpel(:,4))/180*pi;

    EjectTime(i,ipl)=ejecttime/365.25;
    
    SumIList(i,ipl)=abs(sum(di_record_inout));
    StdIList(i,ipl)=(var(di_record_inout))^(1/2);
    RIList(i,ipl)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumeList(i,ipl)=abs(sum(de_record_inout));
    StdeList(i,ipl)=(var(de_record_inout))^(1/2);
    ReList(i,ipl)=sum((de_record_inout-mean(de_record_inout)).^2);
        
end
    
    %% sort IP
    [ePL(:,ipl),ind]=sort(ePL(:,ipl));
    
    aPL(:,ipl)=aPL(ind,ipl);
    aTL(:,ipl)=aTL(ind,ipl);
    eTL(:,ipl)=eTL(ind,ipl);
    IPL(:,ipl)=IPL(ind,ipl);
    ITL(:,ipl)=ITL(ind,ipl);
    
    TimesList(:,ipl)=TimesList(ind,ipl);
    EjectTime(:,ipl)=EjectTime(ind,ipl);
    SumIList(:,ipl)=SumIList(ind,ipl);
    RIList(:,ipl)=RIList(ind,ipl);
    StdIList(:,ipl)=StdIList(ind,ipl);
    SumeList(:,ipl)=SumeList(ind,ipl);
    ReList(:,ipl)=ReList(ind,ipl);
    StdeList(:,ipl)=StdeList(ind,ipl);

    pla0=mean(mean(aPL));
    tpa0=mean(mean(aTL));
    ple0=mean(mean(ePL));
    tpe0=mean(mean(eTL));
    
    e0=1-tpa0/pla0*(1-tpe0);
    e1=1-tpa0/pla0*(1+tpe0);
    
    eval(['e0',num2str(ipl),'=e0;']);
    eval(['e1',num2str(ipl),'=e1;']);

end


end


if ~(exist('NCE','var') && exist('MI','var') && exist('Me','var'))
    
    %% Theo
    NCE=zeros(length(ePL),2);
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

plotx=ePL;

markersize=12;
linewidth=2;
fontsize=12;

xxlim=[0.19 0.41];
yylimNCE=[0.1 1e5];
yylimA=[1e-8 1e0];

figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');

row=4;
col=2;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.40;

subplot(row,col,1);
semilogy(1,1,'w');hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
legend([h1 h2],{'1999CE119&2004UP10','2001FU172&2004UP10'},'fontsize',fontsize/10*8,'location','southeast');

for ipl=1:Npl
        
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    semilogy(plotx(:,ipl),TimesList(:,ipl),[color,'.'],'markersize',markersize);hold on;
    semilogy(plotx(:,ipl),NCE(:,ipl),[color,'.-'],'linewidth',linewidth);

    xlim(xxlim);

end
yylim=yylimNCE;
% ylim(yylim);
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');
hold off;
set(gca,'xticklabel',[]);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
annotation('textbox',[LeftRetainWidth/2.5+Width BottomRetainWidth/3+4*Height 0.03 0.03],'edgecolor','none','string',...,
           '(N)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

       
%% Res time
subplot(row,col,2);
semilogy(0,1,'w');hold all;
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
    xlim(xxlim);

end
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');
set(gca,'xticklabel',[]);
ylabel('$t_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 BottomRetainWidth/3+4*Height 0.03 0.03],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');


subplot(row,col,3);
semilogy(plotx(:,1),SumIList(:,1),'r.','markersize',markersize);hold all;
semilogy(plotx(:,2),SumIList(:,2),'b.','markersize',markersize);
semilogy(plotx(:,1),AI(:,1),'r.-','linewidth',linewidth);
semilogy(plotx(:,2),AI(:,2),'b.-','linewidth',linewidth);
semilogy(plotx(:,1),AISim(:,1),'m--','linewidth',linewidth);
semilogy(plotx(:,2),AISim(:,2),'m--','linewidth',linewidth);

xlim(xxlim);
yylim=yylimA;
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');

hold off;
% set(gca,'xticklabel',[]);
xlabel('$e_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_I|$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth/2.5+Width BottomRetainWidth/3+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(row,col,4);
semilogy(plotx(:,1),SumeList(:,1),'r.','markersize',markersize);hold all;
semilogy(plotx(:,2),SumeList(:,2),'b.','markersize',markersize);
semilogy(plotx(:,1),Ae(:,1),'r.-','linewidth',linewidth);
semilogy(plotx(:,2),Ae(:,2),'b.-','linewidth',linewidth);
semilogy(plotx(:,1),AeSim(:,1),'m--','linewidth',linewidth);
semilogy(plotx(:,2),AeSim(:,2),'m--','linewidth',linewidth);

xlim(xxlim);
yylim=yylimA;
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');
hold off;

% set(gca,'xticklabel',[]);
xlabel('$e_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_e|$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.4+2*Width BottomRetainWidth/3+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
       