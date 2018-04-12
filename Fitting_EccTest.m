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

for ipl=1:Npl
        
    pla=zeros(length(EccList),1);
    tpa=zeros(length(EccList),1);
    ple=zeros(length(EccList),1);
    tpe=zeros(length(EccList),1);

    TimesList=zeros(length(EccList),1);
    SumList=zeros(length(EccList),2);
    StdList=zeros(length(EccList),2);
    RList=zeros(length(EccList),2);
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
for i=1:length(EccList)
    Ecc=EccList(i);
    disp(Ecc);

    fname=[PluName,'_',num2str(sprintf('%.3f',Ecc)),'Ecc'];
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
    pla(i)=mean(plel(:,2));
    ple(i)=mean(plel(:,3));
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
    tpa(i)=mean(tpel(:,2));
    tpe(i)=mean(tpel(:,3));

    TimesList(i)=length(di_record_inout);
    
    SumList(i,1)=abs(sum(di_record_inout));
    StdList(i,1)=(var(di_record_inout))^(1/2);
    RList(i,1)=sum((di_record_inout-mean(di_record_inout)).^2);
    
    SumList(i,2)=abs(sum(de_record_inout));
    StdList(i,2)=(var(de_record_inout))^(1/2);
    RList(i,2)=sum((de_record_inout-mean(de_record_inout)).^2);
        
end
    
    %% sort plotx %% mean value disorder    
    eval(['TimesList',num2str(ipl),'=TimesList;']);
    eval(['SumList',num2str(ipl),'=SumList;']);
    eval(['StdList',num2str(ipl),'=StdList;']);
    eval(['RList',num2str(ipl),'=RList;']);
    
    pla0=mean(pla);
    tpa0=mean(tpa);
    ple0=mean(ple);
    tpe0=mean(tpe);
    
    e0=1-tpa0/pla0*(1-tpe0);
    e1=1-tpa0/pla0*(1+tpe0);
    
    eval(['e0',num2str(ipl),'=e0;']);
    eval(['e1',num2str(ipl),'=e1;']);
    


end
% Last=1;
% EccList(1:Last)=[];
% TimesList(1:Last)=[];
% SumList(1:Last)=[];
% StdList(1:Last)=[];
% RList(1:Last)=[];

% 
% [P,H]=polyfit(log(plotx),log(TimesList),1);  
% R=corrcoef(log(plotx),log(TimesList));
% Timesfit=exp(polyval(P,log(plotx)));
% 
% [P2,H2]=polyfit(log(plotx),log(SumList),1);  
% R2=corrcoef(log(plotx),log(SumList));
% Sumfit=exp(polyval(P2,log(plotx)));
% 
% [P3,H3]=polyfit(log(plotx),log(StdList),1);  
% R3=corrcoef(log(plotx),log(StdList));
% Stdfit=exp(polyval(P3,log(plotx)));
% 


end


plotx=EccList;


% [PRI,HRI]=polyfit(plotx,log(RList(:,1)),1);  
% RRI=corrcoef(plotx,log(RList(:,1)));
% RfitI=exp(polyval(PRI,plotx));


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
plot(1,1,'w');hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'b-');
legend([h1 h2],{'1999CE119&2004UP10','2001FU172&2004UP10'},'fontsize',fontsize/10*8,'location','north');

for ipl=1:Npl
        
    TimesList=eval(['TimesList',num2str(ipl)]);
    e0=eval(['e0',num2str(ipl)]);
    e1=eval(['e1',num2str(ipl)]);
    switch ipl
        case 1 
            color='r';
        case 2
            color='b';
    end
    plot(plotx,TimesList,[color,'.-'],'markersize',markersize);hold on;
    %plot(plotx,Timesfit,'r.-');
    xlim([0.2 0.4]);

%     plot([e0 e0],yylim,'k--','linewidth',2);
%     plot([e1 e1],yylim,'k--','linewidth',2);
end
yylim=[0 6000];
ylim(yylim);
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');
hold off;
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
annotation('textbox',[LeftRetainWidth/2.5+Width BottomRetainWidth/3+4*Height 0.03 0.03],'edgecolor','none','string',...,
           '(A)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

       
subplot(row,col,3);
semilogy(plotx,SumList1(:,1),'r.-','markersize',markersize);hold all;
semilogy(plotx,SumList2(:,1),'b.-','markersize',markersize);
yylim=[1e-8 1e0];
ylim(yylim);
set(gca,'ytick',power(10,-7:-1));
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');

% [AX,H1,H2]=plotyy(plotx,SumList1(:,1),plotx,SumList2(:,1));hold all;
% set(H1,'color','r');
% set(H2,'color','b');
% set(H1,'LineStyle','-');
% set(H2,'LineStyle','-');
% set(H1,'Linewidth',1);
% set(H2,'Linewidth',1);
% set(H1,'markersize',markersize);
% set(H2,'markersize',markersize);
% set(H1,'marker','.');
% set(H2,'marker','.');
% set(AX(1),'ylim',[0 0.22]);
% set(AX(2),'ylim',[0 0.011]);
% set(AX(1),'ytick',0.02:0.02:0.18);
% set(AX(2),'ytick',0.001:0.001:0.009);
% set(AX(1),'XColor','k','YColor','k');
% set(AX(2),'XColor','k','YColor','k');
% yylim=get(AX(1),'ylim');
% plot([e01 e01],yylim,'k--','linewidth',2);
% plot([e11 e11],yylim,'k--','linewidth',2);
% yylim=get(AX(2),'ylim');
% plot([e02 e02],yylim,'k--','linewidth',2);
% plot([e12 e12],yylim,'k--','linewidth',2);
hold off;
% 
% Labels = str2double(get(AX(2),'YTickLabel'));
% scale = 1e3;
% set(AX(2),'YTickLabel',Labels*scale,'units','normalized');
% posAxes = get(AX(2),'position');
% textBox = annotation('textbox','linestyle','none','string',['x 10^{' sprintf('%0.0f',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1)+posAxes(3)/10*8.5 posAxes(2)+posAxes(4)/10*6 posAn(3) posAn(4)],'VerticalAlignment','cap');

set(gca,'xticklabel',[]);
ylabel('$|A_I|=\left|\sum \Delta I\right|~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth/2.5+Width BottomRetainWidth/3+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

       
subplot(row,col,5);
semilogy(plotx,RList1(:,1),'r.-','markersize',markersize);hold all;
semilogy(plotx,RList2(:,1),'b.-','markersize',markersize);
yylim=[1e-8 1.0e-1];
ylim(yylim);
set(gca,'ytick',power(10,-8:-2));
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');

% [AX,H1,H2]=plotyy(plotx,RList1(:,1),plotx,RList2(:,1));hold all;
% yylim1=[0 1.2e-2];
% yylim2=[0 6e-5];
% 
% plot([e01 e01],yylim1,'k--','linewidth',2)
% plot([e11 e11],yylim1,'k--','linewidth',2)
% plot([e02 e02],yylim2,'k--','linewidth',2);
% plot([e12 e12],yylim2,'k--','linewidth',2);
% set(H1,'color','r');
% set(H2,'color','b');
% set(H1,'LineStyle','-');
% set(H2,'LineStyle','-');
% set(H1,'Linewidth',1);
% set(H2,'Linewidth',1);
% set(H1,'markersize',markersize);
% set(H2,'markersize',markersize);
% set(H1,'marker','.');
% set(H2,'marker','.');
% set(AX(1),'XColor','k','YColor','k');
% set(AX(2),'XColor','k','YColor','k');
% set(AX(1),'ylim',yylim1);
% set(AX(2),'ylim',yylim2);
% set(AX(1),'ytick',0:2e-3:1e-2);
% set(AX(1),'yticklabel',sprintfc('%g',0:0.2:1));
% set(AX(2),'ytick',0:1e-5:5e-5);
% 
% Labels = str2double(get(AX(2),'YTickLabel'));
% set(AX(2),'YTickLabel',Labels,'units','normalized');
% scale = 1e5;
% posAxes = get(AX(2),'position');
% textBox = annotation('textbox','linestyle','none','string',['x 10^{' sprintf('%0.0f',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1)+posAxes(3)/10*8.5 posAxes(2)+posAxes(4)/10*5 posAn(3) posAn(4)],'VerticalAlignment','cap');
% 
% scale = 1e-2;
% posAxes = get(AX(1),'position');
% textBox = annotation('textbox','linestyle','none','string',['x 10^{-' sprintf('%0.0f',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1)-posAn(3)/10*6 posAxes(2)+posAn(4)/10*8 posAn(3) posAn(4)],'VerticalAlignment','cap');

% index1=find(EccList<=e11,1,'last');
% fitx1=plotx(index1+1:end);
% [PRe1,HRe1]=polyfit(fitx1,log(RList1(index1+1:end,2)),1);  
% RRe1=corrcoef(fitx1,log(RList1(index1+1:end,2)));
% Rfite1=exp(polyval(PRe1,fitx1));
% 
% index2=find(EccList<=e12,1,'last');
% fitx2=plotx(index2+1:end);
% [PRe2,HRe2]=polyfit(fitx2,log(RList2(index2+1:end,2)),1);  
% RRe2=corrcoef(fitx2,log(RList2(index2+1:end,2)));
% Rfite2=exp(polyval(PRe2,fitx2));
% 
% semilogy(fitx1,Rfite1,'r-','markersize',markersize);
% semilogy(fitx2,Rfite2,'b-','markersize',markersize);

%% quadratic
index0=find(TimesList1>1,1,'first');
fitx=plotx(index0:end);
[PRI,HRI]=polyfit(fitx,log(RList1(index0:end,1)),6);  
RRI1=corrcoef(fitx,log(RList1(index0:end,1)));
RfitI=exp(polyval(PRI,fitx));

semilogy(fitx,RfitI,'r-','markersize',markersize);

index0=find(TimesList2>1,1,'first');
fitx=plotx(index0:end);
[PRI,HRI]=polyfit(fitx,log(RList2(index0:end,1)),6);  
RRI2=corrcoef(fitx,log(RList2(index0:end,1)));
RfitI=exp(polyval(PRI,fitx));

semilogy(fitx,RfitI,'b-','markersize',markersize);


hold off;

xlabel('$e_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
annotation('textbox',[LeftRetainWidth/2.5+Width BottomRetainWidth/3+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(C1)','fontsize',fontsize/10*8,'color','k','fontweight','bold');

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(row,col,4);
%[AX,H1,H2]=plotyy(plotx,SumList1(:,2),plotx,SumList2(:,2));hold all;
semilogy(plotx,SumList1(:,2),'r.-','markersize',markersize);hold all;
semilogy(plotx,SumList2(:,2),'b.-','markersize',markersize);
% set(H1,'color','r');
% set(H2,'color','b');
% set(H1,'LineStyle','-');
% set(H2,'LineStyle','-');
% set(H1,'Linewidth',1);
% set(H2,'Linewidth',1);
% set(H1,'markersize',markersize);
% set(H2,'markersize',markersize);
% set(H1,'marker','.');
% set(H2,'marker','.');
% set(AX(1),'ylim',[0 6e-3]);
% set(AX(2),'ylim',[0 1.2e-3]);
% set(AX(1),'ytick',1e-3:1e-3:6e-3);
% set(AX(2),'ytick',2e-4:2e-4:1.2e-3);
% set(AX(1),'XColor','k','YColor','k');
% set(AX(2),'XColor','k','YColor','k');
% yylim=get(AX(1),'ylim');
% plot([e01 e01],yylim,'k--','linewidth',2);
% plot([e11 e11],yylim,'k--','linewidth',2);
% yylim=get(AX(2),'ylim');
% plot([e02 e02],yylim,'k--','linewidth',2);
% plot([e12 e12],yylim,'k--','linewidth',2);
yylim=[1e-6 1e-1];
ylim(yylim);
set(gca,'ytick',power(10,-5:-1));
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');
hold off;

set(gca,'xticklabel',[]);
ylabel('$|A_e|=\left|\sum \Delta e\right|$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.4+2*Width BottomRetainWidth/3+3*Height 0.03 0.03],'edgecolor','none','string',...,
           '(B2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
       
subplot(row,col,6);
semilogy(plotx,RList1(:,2),'r.-','markersize',markersize);hold all;
semilogy(plotx,RList2(:,2),'b.-','markersize',markersize);

% [AX,H1,H2]=plotyy(plotx,RList1(:,2),plotx,RList2(:,2));hold all;
yylim=[1e-8 1e-4];
% % 
ylim(yylim);
% plot([e01 e01],yylim,'k--','linewidth',2)
% plot([e11 e11],yylim,'k--','linewidth',2)
% plot([e02 e02],yylim,'k--','linewidth',2);
% plot([e12 e12],yylim,'k--','linewidth',2);
patch([e11 e11 e01 e01],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.5,'edgecolor','none');  %???0.5
set(gca,'ytick',power(10,-8:-5));

% set(H1,'color','r');
% set(H2,'color','b');
% set(H1,'LineStyle','-');
% set(H2,'LineStyle','-');
% set(H1,'Linewidth',1);
% set(H2,'Linewidth',1);
% set(H1,'markersize',markersize);
% set(H2,'markersize',markersize);
% set(H1,'marker','.');
% set(H2,'marker','.');
% set(AX(2),'XColor','k','YColor','k');
% set(AX(1),'ylim',yylim1);
% set(AX(2),'ylim',yylim2);
% set(AX(1),'ytick',0:3e-6:1.5e-5);
% set(AX(2),'ytick',0:1e-7:5e-7);
% 
% Labels = str2double(get(AX(2),'YTickLabel'));
% scale = 1e7;
% set(AX(2),'YTickLabel',Labels,'units','normalized');
% posAxes = get(AX(2),'position');
% textBox = annotation('textbox','linestyle','none','string',['x 10^{' sprintf('%0.0f',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1)+posAxes(3) posAxes(2)+posAxes(4)/10*5 posAn(3) posAn(4)],'VerticalAlignment','cap');
% 
% Labels = str2double(get(AX(1),'YTickLabel'));
% scale = 1e5;
% set(AX(1),'YTickLabel',Labels,'units','normalized');
% posAxes = get(AX(1),'position');
% textBox = annotation('textbox','linestyle','none','string',['x 10^{' sprintf('%0.0f',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1)-posAn(3)/10 posAxes(2)+posAxes(4)/10*5 posAn(3) posAn(4)],'VerticalAlignment','cap');
% 
% set(AX(1),'XColor','k','YColor','k');

% %% right
% index1=find(EccList<=roundn(e01,-3),1,'last');
% fitx1=plotx1(index1:end);
% [PRe1,HRe1]=polyfit(fitx1,log(RList1(index1:end,2)),1);  
% RRe1=corrcoef(fitx1,log(RList1(index1:end,2)));
% Rfite1=exp(polyval(PRe1,fitx1));
% 
% index2=find(EccList<=roundn(e02,-3),1,'last');
% fitx2=plotx2(index2:end);
% [PRe2,HRe2]=polyfit(fitx2,log(RList2(index2:end,2)),1);  
% RRe2=corrcoef(fitx2,log(RList2(index2:end,2)));
% Rfite2=exp(polyval(PRe2,fitx2));
% 
% %%% medium
% 
% index1_m=find(EccList>=roundn(e11,-3),1,'first');
% fitx1_m=plotx1(index1_m:index1);
% [PRe1_m,HRe1_m]=polyfit(fitx1_m,log(RList1(index1_m:index1,2)),1);  
% RRe1_m=corrcoef(fitx1_m,log(RList1(index1_m:index1,2)));
% Rfite1_m=exp(polyval(PRe1_m,fitx1_m));
% 
% index2_m=find(EccList>=roundn(e12,-3),1,'first');
% fitx2_m=plotx2(index2_m:index2);
% [PRe2_m,HRe2_m]=polyfit(fitx2_m,log(RList2(index2_m:index2,2)),1);  
% RRe2_m=corrcoef(fitx2_m,log(RList2(index2_m:index2,2)));
% Rfite2_m=exp(polyval(PRe2_m,fitx2_m));
% 
% %%% left
% index0_1=find(TimesList1>1,1,'first');
% fitx1_left=plotx1(index0_1:index1_m);
% [PRe1_left,HRe1_left]=polyfit(fitx1_left,log(RList1(index0_1:index1_m,2)),1);  
% RRe1_left=corrcoef(fitx1_left,log(RList1(index0_1:index1_m,2)));
% Rfite1_left=exp(polyval(PRe1_left,fitx1_left));
% 
% index0_2=find(TimesList2>1,1,'first');
% fitx2_left=plotx2(index0_2:index2_m);
% [PRe2_left,HRe2_left]=polyfit(fitx2_left,log(RList2(index0_2:index2_m,2)),1);  
% RRe2_left=corrcoef(fitx2_left,log(RList2(index0_2:index2_m,2)));
% Rfite2_left=exp(polyval(PRe2_left,fitx2_left));
% 
% semilogy(fitx1,Rfite1,'r-','markersize',markersize);
% semilogy(fitx2,Rfite2,'b-','markersize',markersize);
% 
% semilogy(fitx1_m,Rfite1_m,'r-','markersize',markersize);
% semilogy(fitx2_m,Rfite2_m,'b-','markersize',markersize);
% 
% semilogy(fitx1_left,Rfite1_left,'r-','markersize',markersize);
% semilogy(fitx2_left,Rfite2_left,'b-','markersize',markersize);

%% quadratic
index0=find(TimesList1>1,1,'first');
fitx=plotx(index0:end);
[PRe,HRe]=polyfit(fitx,log(RList1(index0:end,2)),6);  
RRe1=corrcoef(fitx,log(RList1(index0:end,2)));
Rfite=exp(polyval(PRe,fitx));

semilogy(fitx,Rfite,'r-','markersize',markersize);

index0=find(TimesList2>1,1,'first');
fitx=plotx(index0:end);
[PRe,HRe]=polyfit(fitx,log(RList2(index0:end,2)),6);  
RRe2=corrcoef(fitx,log(RList2(index0:end,2)));
Rfite=exp(polyval(PRe,fitx));

semilogy(fitx,Rfite,'b-','markersize',markersize);

% index0_2=find(TimesList2>1,1,'first');
% fitx2_left=plotx(index0_2:index2_m);
% [PRe2_left,HRe2_left]=polyfit(fitx2_left,log(RList2(index0_2:index2_m,2)),1);  
% RRe2_left=corrcoef(fitx2_left,log(RList2(index0_2:index2_m,2)));
% Rfite2_left=exp(polyval(PRe2_left,fitx2_left));


hold off;

xlabel('$e_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2=\sum {(\Delta e-\overline{\Delta e})}^{2}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
annotation('textbox',[LeftRetainWidth*1.4+2*Width BottomRetainWidth/3+2*Height 0.03 0.03],'edgecolor','none','string',...,
           '(C2)','fontsize',fontsize/10*8,'color','k','fontweight','bold');
