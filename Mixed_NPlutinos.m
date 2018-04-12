ffname='RealPlutinos';
fname='1999CE119_1Gyr';

% tpel=load(['~/Documents/ServerMount/LAB Backup/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt']);
% plel=load(['~/Documents/ServerMount/LAB Backup/LAB/CE_realp/',ffname,'/',fname,'/plel.txt']);

% tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt']);
% plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/plel.txt']);

tpel=load(['~/Documents/swiftdata/Trojan/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt']);
plel=load(['~/Documents/swiftdata/Trojan/LAB/CE_realp/',ffname,'/',fname,'/plel.txt']);

temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
if isempty(temp)
    ejectNo=size(tpel,1);
else
    ejectNo=temp;
end
clear temp;

tpel=tpel(1:ejectNo-1,:);
plel=plel(1:ejectNo-1,:);

aP=mean(plel(:,2));
eP=mean(plel(:,3));
IP=mean(plel(:,4))/180*pi;
aT=mean(tpel(:,2));
eT=mean(tpel(:,3));
IT=mean(tpel(:,4))/180*pi;
mP=6.56e-9;

Ru=Fun_Ru(aP,aT,eP,eT,IP,IT);
disp(Ru)

% gm=3.5;
% 
% mP=6.56e-9;
% mN=5.17e-5;
% VSN=(1/aT+mN/aT)^(1/2);
% VS=(1/aT)^(1/2);
% 
% %% determine Ru 
% %[x,y]=solve('x^2/aP^2+y^2/aP^2/(1-eP^2)=1','((x-aP*eP)*cos(IP-IT))^2+y^2=aT^2','x','y');
% [x,y]=solve('(x/cos(IP)+aP*eP)^2/aP^2+y^2/aP^2/(1-eP^2)=1','(x/cos(IT))^2+y^2=aT^2','x','y');
% yy=eval(y);
% for i=1:length(yy)
%     if isreal(yy(i)) && yy(i) > 0
%         WW=yy(i)*2;
%     end
% end
% 
% HH=(abs(aP*(1-eP)*sin(IP)+aT*sin(IT))+abs(aP*(1-eP)*sin(IP)-aT*sin(IT)))/2;
% LL=abs(aP*(1-eP)*cos(IP)-aT*cos(IT));
% RR=(HH*LL*WW)^(1/3);
% TP=1/(1/aP^3)^(1/2);
% Ru=RR/(1e9/2/TP)^(1/2);
% disp('Ru:');
% disp(Ru)
% Vcomp=(mN/aT)^(1/2);
% FRk=@(mu)16/15*mu/(((2/aT-1/aP)^(1/2)-(1/aT)^(1/2))*Vcomp);
FRk=@(nuP)Fun_Rk(aP,aT,nuP);

FRd=@(Np)Ru./Np.^(1/2);

RP=1189;%km
Rs=RP/1.5e8;
FRs=@(nuP)(nuP/1).^(1/3)*Rs;

FRrequired=@(Np)((Ru/FRk(1))^2./Np).^(1/6)*Rs;

Np=exp(log(1):0.01:log(1e4));

figure;
set(gcf,'Position',[400,100,700,500],'color','w');
fontsize=15;

BottomRetainWidth=0.15;
LeftRetainWidth=0.15;
Height=0.7;
Width=0.7;
set(gca,'position',[LeftRetainWidth BottomRetainWidth Width Height]);
annotation('textbox',[(Width+0.06)/2 Height+0.14 0.05 0.05],'edgecolor','none','string',...,
    '1999CE119&2004UP10','fontweight','bold','fontsize',fontsize,'color','k');

annotation('textbox',[LeftRetainWidth+0.01 Height+0.09 0.05 0.05],'edgecolor','none','string',...,
    '$m_{tot}=10\,m_{Pluto}$','Interpreter','latex','fontsize',fontsize,'color','r');

F=@(nu0)((FRk(nu0)-FRs(nu0))*1000);
F1=@(Np)((FRd(Np)-FRs(10/Np))*1000);
Np0=10/fsolve(F,1e-3);
disp('nu0:');
disp(Np0);

h1=loglog(Np,FRk(10./Np),'r-');hold all;
h2=loglog(Np,FRd(Np),'b-');
h3=loglog(Np,FRs(10./Np),'k-','linewidth',2);

%h4=loglog(Np,aP*(10./Np*mP/3).^(1/3),'g-');
% h5=loglog(Np,FRrequired(Np),'k-.','linewidth',2);

yylim=get(gca,'ylim');
xxlim=get(gca,'xlim');
plot([10 10],yylim,'k--');
plot([Np0 Np0],yylim,'k--','linewidth',2);
plot([10 xxlim(2)],[Rs Rs],'--','color',[0.5 0.5 0.5],'linewidth',1);
% fill([Np0 Np0 1e4 1e4],[yylim(1) yylim(2) yylim(2) yylim(1)],'k','facealpha',0.75,'edgecolor','none');
patch([Np(1) Np(1) Np(end) Np(end)],[yylim(1) FRs(10/Np(1)) FRs(10/Np(end)) yylim(1)],'k','facealpha',0.75,'edgecolor','none');
text(Np(floor(length(Np)/3)),FRs(10/Np(floor(length(Np)/3)))/8,'Interior of Plutino','color','w','fontsize',20);

legend([h1 h2 h3 ],{'$R_k$','$R_\varepsilon$','$R_P$'},'fontsize',fontsize,'Interpreter','latex','location','northeast');
xlabel('$N_p$','fontsize',fontsize,'interpreter','latex');
ylabel('$R\rm(AU)$','fontsize',fontsize,'interpreter','latex');
set(gca,'box','off');

lytick=get(gca,'ytick');
ax1=axes('Position',get(gca,'Position'),'Color','none','XAxisLocation','top');
set(ax1,'ytick',[]);
set(ax1,'xscale','log');
set(ax1,'xlim',[Np(1) Np(end)]);
set(ax1,'xticklabel',[])

ax=axes('Position',get(gca,'Position'),'Color','none','YAxisLocation','right');
set(ax,'yscale','log');
lytickRs=lytick/Rs;
set(ax,'ylim',[lytickRs(1) lytickRs(end)]);
%set(ax,'yticklabel',);
set(ax,'xtick',[]);
ylabel('$R/R_{Pluto}$','fontsize',fontsize,'interpreter','latex');

disp('Tau 10:')
disp(FRk(10/10)/FRd(10));
disp('Tau 40:')
disp(FRk(40/40)/FRd(40));