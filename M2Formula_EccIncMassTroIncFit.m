
if exist('TimesListe','var') 
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')

clear;
de_name='de_record_inout';

%% Data collect ------------------------------------------------------------------
%% Ecc ------------------------------------------------------------------
Dir='swiftdata';
DDir='EccTest2';
PluName='1999CE119';

EccList=0.2:0.005:0.35;
EccList=EccList';

ple=zeros(length(EccList),1);
RListe=zeros(length(EccList),1);
TimesListe=zeros(length(EccList),1);

pleplI=zeros(length(EccList),1);
pletpTroI=zeros(length(EccList),1);

for i=1:length(EccList)
    Ecc=EccList(i);
    disp(['pl ecc # ',num2str(Ecc)]);

    fname=[PluName,'_',num2str(sprintf('%.3f',Ecc)),'Ecc'];
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
        
    disp(['  plInc ',num2str(mean(plel(:,4)))]);
    disp(['  tpInc ',num2str(mean(tpel(:,4)))]);
    
    pleplI(i)=mean(plel(:,4));
    pletpTroI(i)=mean(tpel(:,4));
    
    ple(i)=mean(plel(:,3));
    RListe(i)=sum((de_record_inout-mean(de_record_inout)).^2);
    TimesListe(i)=length(de_record_inout);

end
    
%% sort ple %% mean value disorder
[ple,ind]=sort(ple);
RListe=RListe(ind,:);
TimesListe=TimesListe(ind,:);

%%  Inc ------------------------------------------------------------------
Dir='ServerMount';
DDir='IncTest';
PluName='1999CE119';

IncList=0:1:30;
IncList=IncList';

TimesListI=zeros(length(IncList),1);
RListI=zeros(length(IncList),1);
plI=zeros(length(IncList),1);

plIple=zeros(length(IncList),1);
plItpTroI=zeros(length(IncList),1);

for i=1:length(IncList)
    Inc=IncList(i);
    disp(Inc);
    
    fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
    
    plIple(i)=mean(plel(:,3));
    plItpTroI(i)=mean(tpel(:,4));
    
    plI(i)=mean(plel(:,4));
    TimesListI(i)=length(de_record_inout);
    RListI(i)=sum((de_record_inout-mean(de_record_inout)).^2);
end

%% sort pli %% mean value disorder
[plI,ind]=sort(plI);
RListI=RListI(ind,:);
TimesListI=TimesListI(ind,:);

%%  Tro Inc ------------------------------------------------------------------
Dir='ServerMount';
DDir='IncTestTro1';
PluName='2004UP10';

IncList=0:1:30;
IncList=IncList';

TimesListTroI=zeros(length(IncList),1);
RListTroI=zeros(length(IncList),1);
tpTroI=zeros(length(IncList),1);

tpTroIplI=zeros(length(IncList),1);
tpTroIple=zeros(length(IncList),1);

for i=1:length(IncList)
    Inc=IncList(i);
    disp(Inc);
    
    fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
       
    tpTroIplI(i)=mean(plel(:,4));
    tpTroIple(i)=mean(plel(:,3));
    
    tpTroI(i)=mean(tpel(:,4));
    TimesListTroI(i)=length(de_record_inout);
    RListTroI(i)=sum((de_record_inout-mean(de_record_inout)).^2);
end

%% sort pli %% mean value disorder
[tpTroI,ind]=sort(tpTroI);
RListTroI=RListTroI(ind,:);
TimesListTroI=TimesListTroI(ind,:);

%%  Mass ------------------------------------------------------------------
Dir='ServerMount';
PluName='1999CE119';
DDir='MassTest2';

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

if strcmp(PluName,'1999CE119') 
    barrier=5.0;
elseif strcmp(PluName,'1999CE119_2006RJ103') 
    barrier=10.0;
else
    barrier=30.0;
end

iF=find(RatioList2<=barrier,1,'last');
RatioList1=[RatioList1;RatioList2(1:iF)];
RatioList2(1:iF)=[];

TimesList=zeros(length(RatioList),1);
RListm=zeros(length(RatioList),1);

mple=zeros(length(RatioList),1);
mplI=zeros(length(RatioList),1);
mtpTroI=zeros(length(RatioList),1);

for i=1:length(RatioList)
    disp(RatioList(i));
    
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/plel.txt']);
    
    mple(i)=mean(plel(:,3));
    mplI(i)=mean(plel(:,4));
    mtpTroI(i)=mean(tpel(:,4));
    
    %% Caution: ?????1e9 yr??????? ??M2???????
    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if isempty(temp)
        ejectNo=size(tpel,1);
    else
        ejectNo=temp;
    end
    clear temp;
    ejecttime=tpel(ejectNo,1);
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);
    de_record_inout=de_record_inout(1:ejectCENo);  
    
    TimesList(i)=length(de_record_inout);  
    RListm(i)=sum((de_record_inout-mean(de_record_inout)).^2);
end

end

%% Fit and plot -------------------------------------------------

%% Plot confirm
figure;
markersize=15;
fontsize=12;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
row=4;col=1;

%% Zero set
%% must confirm all the data are based on zero-set case!!!
Dir='ServerMount';
DDir='MassTest2';
PluName='1999CE119';
fname=[PluName,'_',num2str(sprintf('%.4f',1.0)),'MP'];
tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

ple0=mean(plel(:,3));
plI0=mean(plel(:,4));
tpTroI0=mean(tpel(:,4));
Ratio0=1.0;
R0=sum((de_record_inout-mean(de_record_inout)).^2);

%% Ecc
subplot(row,col,1);
% yylim=[1e-8 1e-4];
% ylim(yylim);

%% quadratic fit
index0=find(TimesListe>1,1,'first');
fitx=ple(index0:end)/ple0;
fity=RListe(index0:end)/R0;

[Pe,He]=polyfit(fitx,log(fity),6); 
Re=corrcoef(fitx,log(fity));
Rfite=exp(polyval(Pe,fitx));

semilogy(fitx,fity,'k.-','markersize',markersize);hold all;
semilogy(fitx,Rfite,'r-','markersize',markersize);
%ylim([1e-8 1e-4]);
yylim=get(gca,'ylim');
semilogy(ple(1:index0-1)/ple0,[yylim(1)*2 yylim(1)*2],'k.:','markersize',markersize,'linewidth',2);
semilogy([ple(index0-1)/ple0 fitx(1)],[yylim(1)*2 fity(1)],'k.:','markersize',markersize,'linewidth',2);
hold off;

axp=get(gca,'position');
xxlim=get(gca,'xlim');
h=annotation('textarrow',[axp(1)+(ple(index0-1)/ple0-xxlim(1))/(xxlim(2)-xxlim(1))*1.5 ...,
    axp(1)+(ple(index0-1)/ple0-xxlim(1))/(xxlim(2)-xxlim(1))],...,
    [axp(2)+0.017 ...,
    axp(2)+0.017],'string','$M_2=0$',...,
    'Interpreter','latex','textcolor','r','fontsize',fontsize,'color','r',...,
    'headlength',5,'headwidth',5,'headstyle','plain',...,
    'Horizontalalignment','center','verticalalignment','middle');
%(log(yylim(1)*2)-log(yylim(1)))/(log(yylim(2))-log(yylim(1)))]
xlabel('$e_P/e_P^0$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_{2,e}/M_{2,e}^{\,0}$','fontsize',fontsize,'Interpreter','latex');
%xlim([0.2 0.4]);

%set(text(xxlim(1)+(xxlim(2)-xxlim(1))/4*1,yylim(2)/2,'$y=\alpha_0+\alpha_1x+\alpha_2x^2+\alpha_3x^3+\alpha_4x^4+\alpha_5x^5+\alpha_6x^6$'),'Interpreter','latex','fontsize',fontsize,'color','red');
annotation('textbox',[axp(1)+axp(3)/10*3 axp(2)+axp(4)/10*7 0.1 0.04],'string','$y=\alpha_0+\alpha_1x+\alpha_2x^2+\alpha_3x^3+\alpha_4x^4+\alpha_5x^5+\alpha_6x^6$','Interpreter','latex','fontsize',fontsize,'color','red','edgecolor','none');

axes('position',[axp(1)+axp(3)/3 axp(2)+axp(4)/5 axp(3)/4 axp(4)/3.5]);
plot(pletpTroI,pleplI,'k.');hold on;
plot(tpTroI0,plI0,'r.','markersize',20);
xlabel('$I_T~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
ylabel('$I_P~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
hold off;


%% inc 
subplot(row,col,2);

%% Fit
fitx=sind(plI)/sind(plI0);
fity=RListI/R0;
[PI,HI]=polyfit(fitx,log(fity),1);  
RI=corrcoef(fitx,log(fity));
RfitI=exp(polyval(PI,fitx));

semilogy(fitx,fity,'k.-','markersize',markersize);hold all;
semilogy(fitx,RfitI,'r-','markersize',markersize);
hold off;

xlabel('$I_P/I_P^0$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_{2,e}/M_{2,e}^{\,0}$','fontsize',fontsize,'Interpreter','latex');
%xlim([0 0.6]);
%ylim([1e-8 1e-4]);
yylim=get(gca,'ylim');
xxlim=get(gca,'xlim');
% set(text(xxlim(1)+(xxlim(2)-xxlim(1))/10*7,yylim(2)/2,'$y=\alpha_0+\alpha_1x$'),'Interpreter','latex','fontsize',fontsize,'color','red');
axp=get(gca,'position');
annotation('textbox',[axp(1)+axp(3)/10*3 axp(2)+axp(4)/10*7 0.1 0.04],'string','$y=\alpha_0+\alpha_1x$','Interpreter','latex','fontsize',fontsize,'color','red','edgecolor','none');

axes('position',[axp(1)+axp(3)/3*2 axp(2)+axp(4)/5*3 axp(3)/4 axp(4)/3.5]);
plot(plIple,plItpTroI,'k.');hold on;
plot(ple0,tpTroI0,'r.','markersize',20);
xlabel('$e_P$','fontsize',fontsize/3*2,'Interpreter','latex');
ylabel('$I_T~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
hold off;


%% Tro inc 
subplot(row,col,3);

%% Fit
fitx=sind(tpTroI)/sind(tpTroI0);
fity=RListTroI/R0;
[PTroI,HTroI]=polyfit(fitx,log(fity),1);  
RTroI=corrcoef(fitx,log(fity));
RfitTroI=exp(polyval(PTroI,fitx));

semilogy(fitx,fity,'k.-','markersize',markersize);hold all;
semilogy(fitx,RfitTroI,'r-','markersize',markersize);
hold off;

xlabel('$I_T/I_T^0$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_{2,e}/M_{2,e}^{\,0}$','fontsize',fontsize,'Interpreter','latex');
% xlim([0 0.6]);
% ylim([1e-8 1e-4]);
yylim=get(gca,'ylim');
xxlim=get(gca,'xlim');
% set(text(xxlim(1)+(xxlim(2)-xxlim(1))/10*7,yylim(2)/2,'$y=\alpha_0+\alpha_1x$'),'Interpreter','latex','fontsize',fontsize,'color','red');
axp=get(gca,'position');
annotation('textbox',[axp(1)+axp(3)/10*3 axp(2)+axp(4)/10*7 0.1 0.04],'string','$y=\alpha_0+\alpha_1x$','Interpreter','latex','fontsize',fontsize,'color','red','edgecolor','none');

axes('position',[axp(1)+axp(3)/3*2 axp(2)+axp(4)/5*3 axp(3)/4 axp(4)/3.5]);
plot(tpTroIple,tpTroIplI,'k.');hold on;
plot(ple0,plI0,'r.','markersize',20);
xlabel('$e_P$','fontsize',fontsize/3*2,'Interpreter','latex');
ylabel('$I_P~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
hold off;


%% Mass
subplot(row,col,4);

%% Fit
fitx=RatioList1/Ratio0;
fity=RListm(1:length(fitx))/R0;
[Pm,Hm]=polyfit(log(fitx),log(fity),1);  
Rm=corrcoef(log(fitx),log(fity));
Rfitm=exp(polyval(Pm,log(fitx)));

loglog(fitx,fity,'k.-','markersize',markersize);hold all;
loglog(fitx,Rfitm,'r-','markersize',markersize);
hold off;

xlabel('$m/m^0$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_{2,e}/M_{2,e}^{\,0}$','fontsize',fontsize,'Interpreter','latex');
%xlim([9e-4 1e1]);
yylim=get(gca,'ylim');
xxlim=get(gca,'xlim');
% set(text(xxlim(1)+(xxlim(2)-xxlim(1))/10*7,yylim(2)/2,'$y=\alpha_0+\alpha_1x$'),'Interpreter','latex','fontsize',fontsize,'color','red');
axp=get(gca,'position');
annotation('textbox',[axp(1)+axp(3)/10*3 axp(2)+axp(4)/10*7 0.1 0.04],'string','$y=\alpha_0+\alpha_1x$','Interpreter','latex','fontsize',fontsize,'color','red','edgecolor','none');

axes('position',[axp(1)+axp(3)/3*2 axp(2)+axp(4)/5 axp(3)/8 axp(4)/3.5]);
plot(mple,mplI,'k.');hold on;
plot(ple0,plI0,'r.','markersize',20);
xlabel('$e_P$','fontsize',fontsize/3*2,'Interpreter','latex');
ylabel('$I_P~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
hold off;

axes('position',[axp(1)+axp(3)/3*2+axp(3)/8 axp(2)+axp(4)/5 axp(3)/8 axp(4)/3.5]);
plot(mtpTroI,mplI,'k.');hold on;
plot(tpTroI0,plI0,'r.','markersize',20);
set(gca,'ytick',[]);
set(gca,'xtick',[20 40]);
xlabel('$I_T~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
%ylabel('$I_P~\rm(DEG)$','fontsize',fontsize/3*2,'Interpreter','latex');
hold off;

