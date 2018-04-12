% aT=30;aP=40;
% IP=24.694/180*pi;eP=0.272;



if ~exist('TimesList','var')
    clear;
    
    PluNameList={'1999CE119','2001FU172','1999CE119_2006RJ103','2001FU172_2006RJ103'};
    DDirList={'MassTest2','MassTest3','MassTest4','MassTest5'};
    
    for i=1:length(DDirList)
        PluName=PluNameList{i};
        DDir=DDirList{i};

        [aP,aT,eP,eT,IP,IT,RatioList,TimesList,MinDisList,~,~,~]=...,
            Fun_readMassHeap(PluName,DDir);
        eval(['aP',int2str(i),'=aP;']);
        eval(['aT',int2str(i),'=aT;']);
        eval(['eP',int2str(i),'=eP;']);
        eval(['eT',int2str(i),'=eT;']);
        eval(['IP',int2str(i),'=IP;']);
        eval(['IT',int2str(i),'=IT;']);
        eval(['RatioList',int2str(i),'=RatioList;']);
        eval(['TimesList',int2str(i),'=TimesList;']);
        eval(['MinDisList',int2str(i),'=MinDisList;']);
    end
end

%ig=1;

for ig=1:4
    varNameList={'aP','aT','eP','eT','IP','IT',...,
        'RatioList','TimesList','MinDisList'};
    for j=1:length(varNameList)
        varName=varNameList{j};
        eval([varName,'=',varName,int2str(ig),';']);
    end

PluName=PluNameList{ig};
DDir=DDirList{ig};

Ru=Fun_Ru(aP,aT,eP,eT,IP,IT);
disp(Ru)

F=@(nu)((Fun_Rk(aP,aT,nu)-Ru)*1000);

nu0=fsolve(F,10);
disp('nu0:');
disp(nu0);

%% upward NCE
gm=3.5;
mP=6.56e-9;
mN=5.17e-5;
NCEU=@(nu)(gm*aP/Ru)^2*(mP/3)^(2/3)*nu.^(2/3);
%% downward NCE
% NCEL=@(nu)15^2/3^(2/3)/16^2*gm^2*((aP/aT*(2*aP/aT-1))^(1/2)-aP/aT)^2*mN*mP^(-4/3)*nu.^(-4/3);
NCEL=@(nu)1/3^(2/3)/2*gm^2*((aP/aT*(2*aP/aT-1))^(1/2)-aP/aT)^2*mN*mP^(-4/3)*nu.^(-4/3);

figure(ig);
fontsize=15;
set(gcf,'Position',[400,100,500,500*0.618*2],'color','w');

BottomRetainWidth=0.1;
LeftRetainWidth=0.12;
Height=0.42;
Width=0.8;

switch PluName
    case '1999CE119'
        annotation('textbox',[Width/2-0.06 0.08+2*Height 0.05 0.05],'edgecolor','none','string',...,
            '1999CE119&2004UP10','fontweight','bold','fontsize',fontsize,'color','k');
    case '2001FU172'
        annotation('textbox',[Width/2-0.06 0.08+2*Height 0.05 0.05],'edgecolor','none','string',...,
            '2001FU172&2004UP10','fontweight','bold','fontsize',fontsize,'color','k');
    case '1999CE119_2006RJ103'
        annotation('textbox',[Width/2-0.06 0.08+2*Height 0.05 0.05],'edgecolor','none','string',...,
            '1999CE119&2006RJ103','fontweight','bold','fontsize',fontsize,'color','k');
    case '2001FU172_2006RJ103'
        annotation('textbox',[Width/2-0.06  0.08+2*Height 0.05 0.05],'edgecolor','none','string',...,
            '2001FU172&2006RJ103','fontweight','bold','fontsize',fontsize,'color','k');
end

subplot(2,1,1);

loglog(RatioList,TimesList,'k.');hold all;
xxlim=get(gca,'xlim');
xU=xxlim(1):0.1:nu0;
h1=loglog(xU,NCEU(xU),'k-','linewidth',2);
xL=nu0:0.1:xxlim(2);
h2=loglog(xL,NCEL(xL),'k-.','linewidth',2);

legend([h1 h2],{'$N_{CE}(\nu_P<\nu_{P,th})$','$N_{CE}(\nu_P>\nu_{P,th})$'},'fontsize',fontsize,'Interpreter','latex','location','northwest');

set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
set(gca,'xticklabel',[]);
yylim=get(gca,'ylim');
loglog([nu0 nu0],yylim,'k--');
set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));
ylabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

annotation('textbox',[Width+LeftRetainWidth/2 0.08+Height 0.05 0.05],'edgecolor','none','string',...,
           '(A)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

subplot(2,1,2);
loglog(RatioList,MinDisList,'k.');hold all;
%h1=loglog(RatioList,FR(RatioList),'r-');
h1=loglog(RatioList,Fun_Rk(aP,aT,RatioList),'r-','linewidth',2);
% h1=loglog(RatioList,FRS(RatioList*mP),'b-');

xxlim=get(gca,'xlim');
% h2=loglog(xxlim,[FRP(mP) FRP(mP)],'b-');
h2=loglog(xxlim,[Ru Ru],'b-','linewidth',2);
yylim=get(gca,'ylim');
% loglog([nu0 nu0],[yylim(1) FRP(mP)],'k--');
loglog([nu0 nu0],yylim,'k--');

xlabel('$\nu_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$D_{min}~\rm(AU)$','fontsize',fontsize,'Interpreter','latex');
legend([h1 h2],{'$R_k$','$R_u$'},'fontsize',fontsize,'Interpreter','latex','location','northwest');
set(gca,'position',[LeftRetainWidth BottomRetainWidth Width Height]);

yylim=get(gca,'ylim');
set(gca,'yTick',power(10,log10(yylim(1)):1:log10(yylim(2))-1));

annotation('textbox',[Width+LeftRetainWidth/2 0.08 0.05 0.05],'edgecolor','none','string',...,
           '(B)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

end