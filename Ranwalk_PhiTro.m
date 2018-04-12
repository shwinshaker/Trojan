%%Inc test
if ~exist('TimesList','var') 

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

DeltaPhiList=zeros(length(List),2);
TimesList=zeros(length(List),2);

SumIList=zeros(length(List),2);
StdIList=zeros(length(List),2);
RIList=zeros(length(List),2);

SumeList=zeros(length(List),2);
StdeList=zeros(length(List),2);
ReList=zeros(length(List),2);

EjectTime=zeros(length(List),2);


for ipl=1:Npl
        
    DDir=eval(['DDir',num2str(ipl)]);
    TroName=eval(['TroName',num2str(ipl)]);
    
    errid=[];

for i=1:length(List)
    
    disp(['#',num2str(i)]);
    fname=[TroName,'_',num2str(i)];
    
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

    temp=find(tpel(:,2)>30.5 | tpel(:,2)<29.5,1,'first');
    if ~isempty(temp)
        errid=[errid;i];
        continue;
    end
    
    EjectTime(i,ipl)=tpel(end,1)/365.25;
    
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
    di_record_inout=di_record_inout/180*pi;
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt']);

    Nepel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/Nepel.txt']);
    phi=mod((tpel(:,5)+tpel(:,6)+tpel(:,7))-(Nepel(:,5)+Nepel(:,6)+Nepel(:,7)),360);
    DeltaPhiList(i,ipl)=(max(phi)-min(phi))/2;
    
    TimesList(i,ipl)=length(di_record_inout);
    
    SumIList(i,ipl)=abs(sum(di_record_inout));
    StdIList(i,ipl)=std(di_record_inout);
    RIList(i,ipl)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumeList(i,ipl)=abs(sum(de_record_inout));
    StdeList(i,ipl)=std(de_record_inout);
    ReList(i,ipl)=sum((de_record_inout-mean(de_record_inout)).^2);
    
end

    %% sort DeltaPhi
    [DeltaPhiList(:,ipl),ind]=sort(DeltaPhiList(:,ipl));
    
    EjectTime(:,ipl)=EjectTime(ind,ipl);
    TimesList(:,ipl)=TimesList(ind,ipl);
    SumIList(:,ipl)=SumIList(ind,ipl);
    RIList(:,ipl)=RIList(ind,ipl);
    StdIList(:,ipl)=StdIList(ind,ipl);
    SumeList(:,ipl)=SumeList(ind,ipl);
    StdeList(:,ipl)=StdeList(ind,ipl);
    ReList(:,ipl)=ReList(ind,ipl);

end

end

func=@(x)sind(x);
plotf='semilogy';

markersize=12;
fontsize=12;

yylim1=[1e2 1e4];
yylim2=[1e-7 1e-1];

%%  Figure -------------------------------------------------------
figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
row=4;
col=2;

isGridOn=1;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.40;

subplot(row,col,1);
eval([plotf,'(1,1,''w'');']);hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
if isGridOn
    grid on;
end
legend([h1 h2],{'1999CE119&2004UP10','1999CE119&2006RJ103'},...,
    'fontsize',fontsize,'location','southwest');

for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList(:,ipl),TimesList(:,ipl),[color,''.''],''markersize'',markersize);']);

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
           '(N)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

subplot(row,col,2);
eval([plotf,'(0,1,''w'');']);hold all;
if isGridOn
    grid on;
end
for ipl=1:Npl
        
    switch ipl
        case 1
            color='r';
        case 2
            color='b';
    end
    semilogy(DeltaPhiList(:,ipl),EjectTime(:,ipl),[color,'.'],'markersize',markersize);hold all;
    yylim=[1e6 1e10];
    ylim(yylim);
    
end
xlim([0 180]);
set(gca,'xtick',0:30:180);
set(gca,'xticklabel',[]);
ylabel('$t_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 BottomRetainWidth+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
subplot(row,col,3);
eval([plotf,'(1,3,''w'');']);hold all;
if isGridOn
    grid on;
end
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList(:,ipl),SumIList(:,ipl),[color,''.''],''markersize'',markersize);']);
    
end
hold off;

xlim([0 180]);
xlabel('$\Delta\phi~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_I|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xtick',0:30:180);
% set(gca,'xticklabel',[]);
if exist('yylim2','var')
    yylim=yylim2;
else
    yylim=get(gca,'ylim');
end
ylim(yylim);

% yytick=get(gca,'yticklabel');
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*0.5+Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
subplot(row,col,4);
eval([plotf,'(1,3,''w'');']);hold all;
if isGridOn
    grid on;
end
for ipl=1:Npl
    
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    eval([plotf,'(DeltaPhiList(:,ipl),SumeList(:,ipl),[color,''.''],''markersize'',markersize);']);
 
end
hold off;

xlim([0 180]);
xlabel('$\Delta\phi~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$|A_e|$','fontsize',fontsize,'Interpreter','latex');
%set(gca,'box','on');
set(gca,'xtick',0:30:180);
% set(gca,'xticklabel',[]);
if exist('yylim2','var')
    yylim=yylim2;
else
    yylim=get(gca,'ylim');
end
ylim(yylim);

set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.5+2*Width BottomRetainWidth+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
