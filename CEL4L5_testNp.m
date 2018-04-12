%%Inc test

Dir='ServerMount'; 



% DDir1='DeltaPhiTestTro1';
% DDir2='DeltaPhiTestTro2';
% DDir3='DeltaPhiTestTro3';
% DDir4='DeltaPhiTestTro4';
% DDir5='DeltaPhiTestTro5';
% DDir6='DeltaPhiTestTro6';
% DDir7='DeltaPhiTestTro7';

% TroName1='2004UP10';
% TroName2='2006RJ103';
% TroName3=TroName1;
% TroName4='2013KY18';
% TroName5='2010TT191';
% TroName6='2005TN53';
% TroName7='2014QO441';
TroNameList={'2004UP10';'2006RJ103';'2004UP10';...,
    '2013KY18';'2010TT191';'2005TN53';'2014QO441';...,
    '2004KV18';'2011WG157'};

Npl=length(TroNameList);
% color1='r';
% color2='b';
% color3='r';
% color4='g';
% color5='m';
% color6='c';
% color7='k';
colorList={'r';'b';'r';'g';'m';'c';'k';'y';[0.5 0.5 1]};

di_name='di_record_inout';
de_name='de_record_inout';

List=1:100;

ReIplList=[];
for ipl=1:Npl
    if ~exist(['TimesList',num2str(ipl)],'var')
        ReIplList=[ReIplList;ipl];
    end
end

if ~isempty(ReIplList)
    len=length(ReIplList);
    IncNpl(Npl-len+1:end)=[];
    EccNpl(Npl-len+1:end)=[];
    TimesRatioNpl(Npl-len+1:end)=[];
end

if ~exist('IncNpl','var') || ~exist('EccNpl','var') || ~exist('TimesRatioNpl','var')
    IncNpl=zeros(Npl,1);
    EccNpl=zeros(Npl,1);
    TimesRatioNpl=zeros(Npl,1);
    ReIplList=1:Npl;
end

for iRe=1:length(ReIplList)
    
    ipl=ReIplList(iRe);
    disp(['Ipl:  ',num2str(ipl)]);
    
    TimesList=zeros(length(List),1);
    DeltaPhiList=zeros(length(List),1);
    PhiList=zeros(length(List),1);
    IncList=zeros(length(List),1);
    EccList=zeros(length(List),1);
        
%     DDir=eval(['DDir',num2str(ipl)]);
    DDir=['DeltaPhiTestTro',num2str(ipl)];  
    %TroName=eval(['TroName',num2str(ipl)]);
    TroName=TroNameList{ipl};
    
    errid=[];

for i=1:length(List)
    
    disp(['#',num2str(i)]);
    fname=[TroName,'_',num2str(i)];
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);
    try
        plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
        Nepel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/Nepel.txt']);
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
    
    time=Nepel(:,1);
    %phi=mod(3*(plel(:,5)+plel(:,6)+plel(:,7))-2*(Nepel(:,5)+Nepel(:,6)+Nepel(:,7))-(plel(:,5)+plel(:,6)),360);
    try
        phi=mod((tpel(:,5)+tpel(:,6)+tpel(:,7))-(Nepel(:,5)+Nepel(:,6)+Nepel(:,7)),360);
    catch
        errid=[errid;i];
        continue;        
    end
    %% remove wandering ones
    if max(phi)-min(phi)>180
        errid=[errid;i];
        continue;
    end
    
%     figure(1);
%     plot(time,phi,'k.');
%     ylim([0,360]);
%     xlim([0,max(time)]);
    
    DeltaPhiList(i)=(max(phi)-min(phi))/2;
    PhiList(i)=mean(phi);
    TimesList(i)=length(CE_record);
    IncList(i)=mean(tpel(:,4));
    EccList(i)=mean(tpel(:,3));

end

    %% remove err
    DeltaPhiList(errid)=[];
    PhiList(errid)=[];
    IncList(errid)=[];
    EccList(errid)=[];
    TimesList(errid)=[];

%     %% sort DeltaPhi
%     [DeltaPhiList,ind]=sort(DeltaPhiList);    
%     TimesList=TimesList(ind,:);
%     PhiList=PhiList(ind,:);
    
%     %% Add nan to separate tadpole and horseshoe
%     Indsep=find(DeltaPhiList(2:end)-DeltaPhiList(1:end-1)>30);
%     DeltaPhiList=[DeltaPhiList(1:Indsep);nan;DeltaPhiList(Indsep+1:end)];
%     TimesList=[TimesList(1:Indsep);nan;TimesList(Indsep+1:end)];
%     SumList=[SumList(1:Indsep,:);[nan nan];SumList(Indsep+1:end,:)];
%     RList=[RList(1:Indsep,:);[nan nan];RList(Indsep+1:end,:)];
%     StdList=[StdList(1:Indsep,:);[nan nan];StdList(Indsep+1:end,:)];
    
    % eval(['tpIncMeanData',num2str(ipl),'=mean(tpInc);']);
    eval(['TimesList',num2str(ipl),'=TimesList;']);
    eval(['DeltaPhiList',num2str(ipl),'=DeltaPhiList;']);
    eval(['PhiList',num2str(ipl),'=PhiList;']);
    eval(['EccList',num2str(ipl),'=EccList;']);
    eval(['IncList',num2str(ipl),'=IncList;']);

    IncNpl(ipl)=mean(IncList);
    EccNpl(ipl)=mean(EccList);
    TimesRatioNpl(ipl)=mean(TimesList(PhiList>180))/mean(TimesList(PhiList<180));
    
end

func=@(x)sind(x);
plotf='semilogy';
% plotf='plot';

% plotx=func(IncList);
% tpIncMean1=func(tpIncMeanData1);
% tpIncMean2=func(tpIncMeanData2);

%xtick=[30 20 10 5 3 2 1 0.5 0.3 0.2 0.1];

markersize=15;
fontsize=15;

yylim1=[1e2 1e4];
yylim4=[1e-7 1];
yylim6=[1e-8 1e-3];

%%  Figure -------------------------------------------------------
figure;
% set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
% row=4;
% col=2;
semilogy(0,1,'w+');hold all;

for ipl=1:Npl
    
    DeltaPhiList=eval(['DeltaPhiList',num2str(ipl)]);
    PhiList=eval(['PhiList',num2str(ipl)]);
    TimesList=eval(['TimesList',num2str(ipl)]);
    % tpIncMean=eval(['tpIncMean',num2str(ipl)]);
    %eval(['color=color',num2str(ipl),';']);
    color=colorList{ipl};
%     eval([plotf,'(DeltaPhiList,TimesList(:,1),[color,''.-''],''markersize'',markersize);']);
    eval([plotf,'(PhiList,TimesList(:),''.'',''color'',color,''markersize'',markersize);']);

%     yylim=[1e2 1e4];
%     ylim(yylim);

%     plot([tpIncMean tpIncMean],yylim,'k--','linewidth',2);
    
end
hold off;

xlim([0 360]);
set(gca,'xtick',0:30:360);
xlabel('$\Delta \phi$','Interpreter','Latex','fontsize',fontsize);
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

figure;
set(gcf,'Position',[400,100,800,400],'color','w');

subplot(1,2,1);
plot(IncNpl,TimesRatioNpl,'k.','markersize',markersize);hold on;
axis([0 30 0.5 1.5]);
xxlim=get(gca,'xlim');
plot(xxlim,[1 1],'k--');
xlabel('$Inc.$','Interpreter','Latex','fontsize',fontsize);
ylabel('$N_{L5}/N_{L4}$','Interpreter','Latex','fontsize',fontsize);

subplot(1,2,2);
plot(EccNpl,TimesRatioNpl,'k.','markersize',markersize);hold on;
axis([0 0.5 0.5 1.5]);
xxlim=get(gca,'xlim');
plot(xxlim,[1 1],'k--');
xlabel('$Ecc.$','Interpreter','Latex','fontsize',fontsize);
ylabel('$N_{L5}/N_{L4}$','Interpreter','Latex','fontsize',fontsize);

% figure;
% set(gcf,'Position',[400,100,800,400],'color','w');
% 
% subplot(1,2,1);
% semilogy(0,1,'w+');hold all;
% for ipl=1:Npl
%     
%     EccList=eval(['EccList',num2str(ipl)]);
%     IncList=eval(['IncList',num2str(ipl)]);
%     PhiList=eval(['PhiList',num2str(ipl)]);
%     TimesList=eval(['TimesList',num2str(ipl)]);
%     color=eval(['color',num2str(ipl)]);
%     
%     plot(IncList,mean(TimesList(PhiList>180))/mean(TimesList(PhiList<180)),'.','color',color,'markersize',markersize);hold on;
%     
% end
% %axis([0 30 0.5 1.5]);
% xxlim=get(gca,'xlim');
% plot(xxlim,[1 1],'k--');
% xlabel('$Inc.$','Interpreter','Latex','fontsize',fontsize);
% ylabel('$N_{L5}/N_{L4}$','Interpreter','Latex','fontsize',fontsize);


