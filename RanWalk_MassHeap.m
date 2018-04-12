%%MassTest

if ~exist('TimesList','var')
    clear;
    
    PluNameList={'1999CE119','2001FU172','1999CE119_2006RJ103','2001FU172_2006RJ103'};
    DDirList={'MassTest2','MassTest3','MassTest4','MassTest5'};
    
    for i=1:length(DDirList)
        PluName=PluNameList{i};
        DDir=DDirList{i};

        [aP,aT,eP,eT,IP,IT,RatioList,TimesList,~,SumList,RList,MaxCumsumList,EjectTime]=...,
            Fun_readMassHeap(PluName,DDir);
        varNameList={'aP','aT','eP','eT','IP','IT',...,
            'RatioList','TimesList','SumList','RList','MaxCumsumList','EjectTime'};
        for j=1:length(varNameList)
            varName=varNameList{j};
            eval([varName,int2str(i),'=',varName,';']);
        end            
    end
end

%% Choose pair
ig=1; %% 1999CE119
varNameList={'aP','aT','eP','eT','IP','IT',...,
    'RatioList','TimesList','SumList','RList','MaxCumsumList','EjectTime'};
for j=1:length(varNameList)
    varName=varNameList{j};
    eval([varName,'=',varName,int2str(ig),';']);
end
PluName=PluNameList{ig};

%% Determine mass threshold
Ru=Fun_Ru(aP,aT,eP,eT,IP,IT);
% disp(Ru)
F=@(nu)((Fun_Rk(aP,aT,nu)-Ru)*1000);
nu0=fsolve(F,10);
disp('nu0:');
disp(nu0);

iF=find(RatioList<=nu0,1,'last');
RatioListU=RatioList(1:iF);
RatioListL=RatioList(iF+1:end);

%% Theo
%% upward NCE
gm=3.5;
mP=6.56e-9;
mN=5.17e-5;
Fun_NCEU=@(nu,aP)(gm*aP/Ru)^2*(mP/3)^(2/3)*nu.^(2/3);
%% downward NCE
% Fun_NCEL=@(nu,aP,aT)15^2/3^(2/3)/16^2*gm^2*((aP/aT*(2*aP/aT-1))^(1/2)-aP/aT)^2*mN*mP^(-4/3)*nu.^(-4/3);
Fun_NCEL=@(nu,aP,aT)1/3^(2/3)/2*gm^2*((aP/aT*(2*aP/aT-1))^(1/2)-aP/aT)^2*mN*mP^(-4/3)*nu.^(-4/3);

NCEU=Fun_NCEU(RatioListU,aP);
NCEL=Fun_NCEL(RatioListL,aP,aT);
NCE=[NCEU;NCEL];
%% Extend
NCEUL=Fun_NCEU(RatioListL,aP);

[dinc,decc]=Fun_diDstb_theo(100000,1,aP,eP,IP/pi*180,aT,eT,IT/pi*180);
dinc=dinc/180*pi;
MI=(std(dinc)*RatioList.^(2/3)).^2.*NCE;
AI=sqrt(2/pi*MI);
%% Extend
MIUL=(std(dinc)*RatioListL.^(2/3)).^2.*NCEUL;
AIUL=sqrt(2/pi*MIUL);

% [~,de]=Fun_diDstb_theo(100000,1,aP,eP,IP/pi*180,aT,eT,IT/pi*180);
Me=(std(decc)*RatioList.^(2/3)).^2.*NCE;
Ae=sqrt(2/pi*Me);
%% Extend
MeUL=(std(decc)*RatioListL.^(2/3)).^2.*NCEUL;
AeUL=sqrt(2/pi*MeUL);

%% SimTheo
%% Calculate ChiN and ChiM
ChiN=(3.5*aP/Ru)^2*(mP/3)^(2/3);
disp('lnChiN');
disp(log(ChiN));
A=(3-aT/aP)/2;
B=(2-aT/aP)^(1/2);
ChiM=(aT/Ru)^2*(2/(A-B*cos(max(IP,IT))))*mP^2;
disp('lnChiM');
disp(log(ChiM));
MSim=@(nuP)ChiM*nuP.^2.*(-1+log(ChiN)+log(nuP)/3*2);
MSimTheo=MSim(RatioList);
% ASimTheo=sqrt(2/pi*MSimTheo);
ASimTheo=Fun_totEff_sim(RatioList,aP,eP,IP,aT,eT,IT);

%% Plot
linewidth=2;
fontsize=12;

xxlim=[1e-3 1e3];
xxtick=power(10,-3:1:3);
xp=1.5e-3;

ReviseLineStyle='-';
ReviseColor='c';

figure;
BottomRetainWidth=0.05;
LeftRetainWidth=0.1;
Height=0.23;
Width=0.39;

%% lim
if strcmp(PluName,'1999CE119')   
%     yylimAI=[1e-5 1e3];
    yylimAI=[1e-7 1e1];
else
    yylimAI=[1e-7 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimAe=[1e-7 1e1];
else
    yylimAe=[1e-8 1e1];
end

if strcmp(PluName,'1999CE119')   
%     yylimGI=[1e-5 1e3];
    yylimGI=[1e-7 1e1];
else
    yylimGI=[1e-7 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimGe=[1e-7 1e1];
else
    yylimGe=[1e-8 1e1];
end

if strcmp(PluName,'1999CE119')   
%     yylimMI=[1e-10 1e5];
    yylimMI=[1e-13 1e1];
else
    yylimMI=[1e-12 1e3];
end
if strcmp(PluName,'1999CE119')   
    yylimMe=[1e-13 1e1];
else
    yylimMe=[1e-15 1e1];
end

labelx='$\nu_P$';
ReviseNo=iF+1;

set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');

row=4;
col=2;

col1=[col+1 2*col+1 3*col+1]';
col2=col1+1;
col3=col2+1;

subplot(row,col,1);
TimesReviseList=TimesList(ReviseNo:end,1)./EjectTime(ReviseNo:end)*1e9;
loglog(RatioList,TimesList,'k.');hold all;
loglog(RatioListL,TimesReviseList,'c.');
loglog(RatioListL,NCEUL,'r-.','linewidth',linewidth);
%% Theo
% loglog(RatioListU,NCEU,'r-','linewidth',linewidth);
% loglog(RatioListL,NCEL,'r-','linewidth',linewidth);
loglog(RatioList,NCE,'r-','linewidth',linewidth);

%% fit
[PT1,HT1]=polyfit(log(RatioListU),log(TimesList(1:iF,1)),1);  
RT1=corrcoef(log(RatioListU),log(TimesList(1:iF,1)));

[PT2,HT2]=polyfit(log(RatioListL),log(TimesList(iF+1:end,1)),1);  
RT2=corrcoef(log(RatioListL),log(TimesList(iF+1:end,1)));

if strcmp(PluName,'1999CE119')   
    yylim=[10 1e6];
else
    yylim=[1 1e5];
end
ylim(yylim);
xlim(xxlim);
plot([nu0 nu0],yylim,'r--');
hold off;

set(gca,'xTick',xxtick);
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

set(text(xp,5/10*yylim(2),['$$\ln{y} = ',num2str(PT1(1),'%.4f'),'\,\ln{x}+',num2str(PT1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,2/10*yylim(2),['$$R^2 = ',num2str(RT1(1,2)*RT1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(N)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% Res time
subplot(row,col,2);
loglog(RatioList,EjectTime,'k.');hold all;
% loglog(RatioList,MaxTime(:,1),'r.-');
% loglog(RatioList,MaxTime(:,2),'b.-');
% loglog(RatioList,MaxTime(:,3),'g.-');
yylim=[1e5 1e10];
ylim(yylim);
xlim(xxlim);
xxlim=get(gca,'xlim');
plot([nu0 nu0],yylim,'r--');
ylabel('$T_{res}~\mathrm{(yr)}$','fontsize',fontsize,'Interpreter','latex');
set(gca,'xTick',[]);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));

set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(T)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% Ai
subplot(row,col,col1(1));
loglog(RatioList,SumList(:,1),'k.');hold all;
loglog(RatioList,ASimTheo,'m--','linewidth',linewidth);
loglog(RatioList,AI,'r-','linewidth',linewidth);

%% extend
loglog(RatioListL,AIUL,'r-.','linewidth',linewidth);
plot([nu0 nu0],yylimAI,'r--');
hold off;

ylabel('$A_I$','fontsize',fontsize,'Interpreter','latex');
yylim=yylimAI;
ylim(yylim);
xlim(xxlim);

set(gca,'xTick',xxtick);
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

[PSum1,HSum1]=polyfit(log(RatioListU),log(SumList(1:iF,1)),1);  
RSum1=corrcoef(log(RatioListU),log(SumList(1:iF,1)));

[PSum2,HSum2]=polyfit(log(RatioListL),log(SumList(iF+1:end,1)),1);  
RSum2=corrcoef(log(RatioListL),log(SumList(iF+1:end,1)));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSum1(1),'%.4f'),'\,\ln{x}',num2str(PSum1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RSum1(1,2)*RSum1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
%% cumsum i
subplot(row,col,col1(2));
loglog(RatioList,MaxCumsumList(:,1),'k.');hold all;
loglog(RatioList,2*ASimTheo,'m--','linewidth',linewidth);
loglog(RatioList,2*AI,'r-','linewidth',linewidth);
loglog(RatioListL,2*AIUL,'r-.','linewidth',linewidth);

plot([nu0 nu0],yylimGI,'r--');
hold off;

ylabel('$G_I$','fontsize',fontsize,'Interpreter','latex');
yylim=yylimGI;
ylim(yylim);
xlim(xxlim);

set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

[PCumsum1,HCumsum1]=polyfit(log(RatioListU),log(MaxCumsumList(1:iF,1)),1);  
RCumsum1=corrcoef(log(RatioListU),log(MaxCumsumList(1:iF,1)));

[PCumsum2,HCumsum2]=polyfit(log(RatioListL),log(MaxCumsumList(iF+1:end,1)),1);  
RCumsum2=corrcoef(log(RatioListL),log(MaxCumsumList(iF+1:end,1)));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsum1(1),'%.4f'),'\,\ln{x}',num2str(PCumsum1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsum1(1,2)*RCumsum1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(B1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% M2
subplot(row,col,col1(3));
loglog(RatioList,RList(:,1),'k.');hold all;
loglog(RatioList,MSimTheo,'m--','linewidth',linewidth);
loglog(RatioList,MI,'r-','linewidth',linewidth);
loglog(RatioListL,MIUL,'r-.','linewidth',linewidth);

plot([nu0 nu0],yylimMI,'r--');
hold off;

yylim=yylimMI;
ylim(yylim);
xlim(xxlim);
set(gca,'xTick',xxtick);
set(gca,'yTick',power(10,log10(yylim(1)):2:log10(yylim(2))-1));

[PR1,HR1]=polyfit(log(RatioListU),log(RList(1:iF,1)),1);  
RR1=corrcoef(log(RatioListU),log(RList(1:iF,1)));

[PR2,HR2]=polyfit(log(RatioListL),log(RList(iF+1:end,1)),1);  
RR2=corrcoef(log(RatioListL),log(RList(iF+1:end,1)));

set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PR1(1),'%.4f'),'\,\ln{x}',num2str(PR1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,3/1000*yylim(2),['$$R^2 = ',num2str(RR1(1,2)*RR1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');

xlabel(labelx,'fontsize',fontsize,'Interpreter','latex');
% ylabel('$M_I=\sum {(\Delta I-\overline{\Delta I})}^{2}~\mathrm{({DEG}^2)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_I$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[LeftRetainWidth BottomRetainWidth Width Height]);
annotation('textbox',[Width+LeftRetainWidth/2 0.025 0.05 0.05],'edgecolor','none','string',...,
           '(C1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%de%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(row,col,col2(1));
loglog(RatioList,SumList(:,2),'k.');hold all;
loglog(RatioList,ASimTheo,'m--','linewidth',linewidth);
loglog(RatioList,Ae,'r-','linewidth',linewidth);
%% extend
loglog(RatioListL,AeUL,'r-.','linewidth',linewidth);
plot([nu0 nu0],yylimAe,'r--');
hold off;

ylabel('$A_e$','fontsize',fontsize,'Interpreter','latex');

yylim=yylimAe;
ylim(yylim);
xlim(xxlim);

set(gca,'xTick',xxtick);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));

[PSumde1,HSumde1]=polyfit(log(RatioListU),log(SumList(1:iF,2)),1);  
RSumde1=corrcoef(log(RatioListU),log(SumList(1:iF,2)));

[PSumde2,HSumde2]=polyfit(log(RatioListL),log(SumList(iF+1:end,2)),1);  
RSumde2=corrcoef(log(RatioListL),log(SumList(iF+1:end,2)));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PSumde1(1),'%.4f'),'\,\ln{x}',num2str(PSumde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RSumde1(1,2)*RSumde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

       
%%%%%% G
subplot(row,col,col2(2));
loglog(RatioList,MaxCumsumList(:,2),'k.');hold all;
loglog(RatioList,2*ASimTheo,'m--','linewidth',linewidth);
loglog(RatioList,2*Ae,'r-','linewidth',linewidth);
loglog(RatioListL,2*AeUL,'r-.','linewidth',linewidth);

ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');

yylim=yylimGe;
ylim(yylim);
xlim(xxlim);

plot([nu0 nu0],yylim,'r--');
hold off;
set(gca,'xTick',xxtick);

set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))-1));

[PCumsumde1,HCumsumde1]=polyfit(log(RatioListU),log(MaxCumsumList(1:iF,2)),1);  
RCumsumde1=corrcoef(log(RatioListU),log(MaxCumsumList(1:iF,2)));

[PCumsumde2,HCumsumde2]=polyfit(log(RatioListL),log(MaxCumsumList(iF+1:end,2)),1);  
RCumsumde2=corrcoef(log(RatioListL),log(MaxCumsumList(iF+1:end,2)));

set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsumde1(1),'%.4f'),'\,\ln{x}',num2str(PCumsumde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsumde1(1,2)*RCumsumde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(gca,'xticklabel',[]);
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(B2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

%% M
subplot(row,col,col2(3));
loglog(RatioList,RList(:,2),'k.');hold all;
loglog(RatioList,MSimTheo,'m--','linewidth',linewidth);
loglog(RatioList,Me,'r-','linewidth',linewidth);
loglog(RatioListL,MeUL,'r-.','linewidth',linewidth);

yylim=yylimMe;

ylim(yylim);
xlim(xxlim);
plot([nu0 nu0],yylim,'r--');
hold off;
set(gca,'yTick',power(10,log10(yylim(1)):2:log10(yylim(2))-1));

[PRde1,HRde1]=polyfit(log(RatioListU),log(RList(1:iF,2)),1);  
RRde1=corrcoef(log(RatioListU),log(RList(1:iF,2)));

[PRde2,HRde2]=polyfit(log(RatioListL),log(RList(iF+1:end,2)),1);  
RRde2=corrcoef(log(RatioListL),log(RList(iF+1:end,2)));

set(text(xp,8/100*yylim(2),['$$\ln{y} = ',num2str(PRde1(1),'%.4f'),'\,\ln{x}',num2str(PRde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,3/1000*yylim(2),['$$R^2 = ',num2str(RRde1(1,2)*RRde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');

ylabel('$M_e$','fontsize',fontsize,'Interpreter','latex');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth Width Height]);
annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025 0.05 0.05],'edgecolor','none','string',...,
           '(C2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
set(gca,'xTick',xxtick);
xlabel(labelx,'fontsize',fontsize,'Interpreter','latex');

