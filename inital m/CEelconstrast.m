clear;
ffname='RealPlutinosNpl';
ffname2='RanPlutinos';
%fname='2001FU172&2006RJ103';
%fname='1999CE119&2006RJ103';
fname='1999CE119';
%fname='2001FU172';

Dir=[fname,'_1Gyr_40pl'];
Dir2=[fname,'_1Gyr_40pl_ran'];
%Dir2=[fname,'_40plTimes'];
AE_tp=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',Dir,'/AE_record_pl.txt']);
AE_tp_CE=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',Dir2,'/AE_record_pl.txt']);

inclim=[15 40];
figure(1);
fontsize=15;
title1='elements distribution in total life time';
title2='elements distribution only in CE';
set(gcf,'Position',[400,100,600,600],'color','w');
annotation(gcf,'textbox','String',{title1},'FontSize',fontsize,'Position',[0.30 0.87 0.10 0.10],'edgecolor',get(gcf,'color'))
annotation(gcf,'textbox','String',{title2},'FontSize',fontsize,'Position',[0.30 0.40 0.10 0.10],'edgecolor',get(gcf,'color'))

subplot(2,2,3);
plot(AE_tp_CE(:,1),AE_tp_CE(:,3),'k.');
axis square;xlabel('a (AU)','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
%xlim([39 40]);ylim(inclim);

xlim1=get(gca,'xlim');
ylim1=get(gca,'ylim');

subplot(2,2,4);
plot(AE_tp_CE(:,2),AE_tp_CE(:,3),'k.');
axis square;xlabel('e','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
%xlim([0.15 0.4]);ylim(inclim);
% 
xlim2=get(gca,'xlim');
ylim2=get(gca,'ylim');

subplot(2,2,1);
plot(AE_tp(:,1),AE_tp(:,3),'k.');
axis square;xlabel('a (AU)','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
%xlim([39 40]);ylim(inclim);
axis([xlim1 ylim1]);

subplot(2,2,2);
plot(AE_tp(:,2),AE_tp(:,3),'k.');
axis square;xlabel('e','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
%xlim([0.15 0.4]);ylim(inclim);
axis([xlim2 ylim2]);
