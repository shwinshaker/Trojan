clear;
ffname='RealPlutinosNpl';
%ffname2='RanPlutinos';
%fname='2001FU172&2006RJ103';
fname='1999CE119&2006RJ103';
%fname='1999CE119';
%fname='2001FU172';

MountDir='swiftdata';

Dir=[fname,'_1Gyr_40pl'];
%Dir2=[fname,'_1Gyr_40pl_ran'];
%Dir2=[fname,'_40plTimes'];
AE_pl=load(['~/Documents/',MountDir,'/LAB/CE_realp/',ffname,'/',Dir,'/AE_record_pl.txt']);
%AE_tp_CE=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',Dir2,'/AE_record_pl.txt']);
plCENo=load(['~/Documents/',MountDir,'/LAB/CE_realp/',ffname,'/',Dir,'/plCENo.txt']);



fontsize=15;
inclim=[0 10];
% title1='elements distribution in total life time';
% title2='elements distribution only in CE';
% set(gcf,'Position',[400,100,600,600],'color','w');
% annotation(gcf,'textbox','String',{title1},'FontSize',fontsize,'Position',[0.30 0.87 0.10 0.10],'edgecolor',get(gcf,'color'))
% annotation(gcf,'textbox','String',{title2},'FontSize',fontsize,'Position',[0.30 0.40 0.10 0.10],'edgecolor',get(gcf,'color'))

figure(1);
% subplot(1,2,1);
plot(0,0,'w');hold all;
xlabel('a (AU)','fontsize',fontsize,'Interpreter','latex');
ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
for ip=3:42 
    itemp=find(plCENo==ip);
    plot(AE_pl(itemp,1),AE_pl(itemp,3),'.','color',[(ip-3)/39 0 1-(ip-3)/39]);
    axis square;
    xlim([39 40]);ylim(inclim);
    pause(0.5);
end
hold off;


% xlim1=get(gca,'xlim');
% ylim1=get(gca,'ylim');

% subplot(1,2,2);
% plot(AE_tp_CE(:,2),AE_tp_CE(:,3),'k.');
% axis square;xlabel('e','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
% %xlim([0.15 0.4]);ylim(inclim);
% % 
% xlim2=get(gca,'xlim');
% ylim2=get(gca,'ylim');
% 
% subplot(2,2,1);
% plot(AE_tp(:,1),AE_tp(:,3),'k.');
% axis square;xlabel('a (AU)','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
% %xlim([39 40]);ylim(inclim);
% axis([xlim1 ylim1]);
% 
% subplot(2,2,2);
% plot(AE_tp(:,2),AE_tp(:,3),'k.');
% axis square;xlabel('e','fontsize',fontsize,'Interpreter','latex');ylabel('inc (DEG)','fontsize',fontsize,'Interpreter','latex');
% %xlim([0.15 0.4]);ylim(inclim);
% axis([xlim2 ylim2]);
