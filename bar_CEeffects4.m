%% Four diff. particles

clear;

ffname='RealPlutinos';

fname1='1999CE119_1Gyr';
fname2='2001FU172_1Gyr'; 
fname3='1999CE119&2006RJ103_1Gyr'; 
fname4='2001FU172&2006RJ103_1Gyr';

name1='1999CE119 & 2004UP10';
name2='2001FU172 & 2004UP10';
name3='1999CE119 & 2006RJ103';
name4='2001FU172 & 2006RJ103';

tag='di';

% if strcmp(tag,'de')
%     lim=0.0005;
% else
%     lim=0.01;
% end

xxlim=4.0;

fontsize=15;
Row=2;
Column=2;

figure(1);
set(gcf,'Position',[400,100,600,600],'color','w');
%annotation(gcf,'textbox','String',{'record-inout'},'FontSize',12,'Position',[0.45 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))


for isub=1:4
subplot(Row,Column,isub);
fname=eval(['fname',num2str(isub)]);
name=eval(['name',num2str(isub)]);
CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt'));
PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));

hill=sqrt(r2hill_record);
PV_rlt=PV_record_tp-PV_record_pl;
x=(PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2);
y=PV_rlt(:,3);
theta=abs(atan(y./x)/pi*180);
%theta=atan(y./x)/pi*180;
plot(0,0,'w');hold all;axis square;
    for i=1:length(CE_record)
        plot(CE_record(i,2)/hill,eval(strcat(tag,'_record_inout(i)')),'.','color',[((theta(i)/90)) 0 1-((theta(i)/90))]);
    end
hold off;
grid on;
if mod(isub,2)==1
    yylim=0.01;
else
    yylim=0.001;
end
set(text(1.5,3/4*yylim,['$sum=',num2str(sum(eval(strcat(tag,'_record_inout'))),4),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
set(text(1.5,2.4/4*yylim,['$N_{CE}= ',num2str(length(theta)),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
xlim([0 xxlim]);
ylim([-yylim yylim]);
xlabel('$Dis.~\rm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\delta Inc~\rm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
title(name);
clear 'CE_record' ;'r2hill_record' ;strcat(tag,'_record_inout');
end
% 
% subplot(Row,Column,2);
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/CE_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/r2hill_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/',tag,'_record_inout.txt'));
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/PV_record_tp.txt'));
% 
% hill=sqrt(r2hill_record);
% PV_rlt=PV_record_tp-PV_record_pl;
% x=(PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2);
% y=PV_rlt(:,3);
% theta=abs(atan(y./x)/pi*180);
% 
% plot(0,0,'w');hold all;axis square;
% for i=1:length(CE_record)
%     plot(CE_record(i,2)/hill,eval(strcat(tag,'_record_inout(i)')),'.','color',[((theta(i)/90)) 0 1-((theta(i)/90))]);
% end
% hold off;
% grid on;
% set(text(1.5,3/4*lim,['$Sum=',num2str(sum(eval(strcat(tag,'_record_inout'))),4),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% set(text(1.5,2.4/4*lim,['$N_{CE}= ',num2str(length(theta)),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% xlim([0 xxlim]);
% ylim([-lim lim]);
% xlabel('$Dis.~\rm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\delta Inc~\rm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% title(name2);
% clear 'CE_record';'r2hill_record';strcat(tag,'_record_inout');
% 
% subplot(Row,Column,3);
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname3,'/CE_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname3,'/r2hill_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname3,'/',tag,'_record_inout.txt'));
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname3,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname3,'/PV_record_tp.txt'));
% 
% hill=sqrt(r2hill_record);
% PV_rlt=PV_record_tp-PV_record_pl;
% x=(PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2);
% y=PV_rlt(:,3);
% theta=abs(atan(y./x)/pi*180);
% 
% plot(0,0,'w');hold all;axis square;
% for i=1:length(CE_record)
%     plot(CE_record(i,2)/hill,eval(strcat(tag,'_record_inout(i)')),'.','color',[((theta(i)/90)) 0 1-((theta(i)/90))]);
% end
% hold off;
% grid on;
% set(text(1.5,3/4*lim,['$Sum=',num2str(sum(eval(strcat(tag,'_record_inout'))),4),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% set(text(1.5,2.4/4*lim,['$N_{CE}= ',num2str(length(theta)),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% xlim([0 xxlim]);
% ylim([-lim lim]);
% xlabel('$Dis.~\rm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\delta Inc~\rm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% title(name3);
% clear 'CE_record' ;'r2hill_record' ;strcat(tag,'_record_inout');
% 
% subplot(Row,Column,4);
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname4,'/CE_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname4,'/r2hill_record.txt'));
% load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname4,'/',tag,'_record_inout.txt'));
% PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname4,'/PV_record_pl.txt'));
% PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname4,'/PV_record_tp.txt'));
% 
% hill=sqrt(r2hill_record);
% PV_rlt=PV_record_tp-PV_record_pl;
% x=(PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2);
% y=PV_rlt(:,3);
% theta=abs(atan(y./x)/pi*180);
% 
% plot(0,0,'w');hold all;axis square;
% for i=1:length(CE_record)
%     plot(CE_record(i,2)/hill,eval(strcat(tag,'_record_inout(i)')),'.','color',[((theta(i)/90)) 0 1-((theta(i)/90))]);
% end
% hold off;
% grid on;
% set(text(1.5,3/4*lim,['$Sum=',num2str(sum(eval(strcat(tag,'_record_inout'))),4),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% set(text(1.5,2.4/4*lim,['$N_{CE}= ',num2str(length(theta)),'$']),'color','red','fontsize',fontsize,'Interpreter','latex');
% 
% xlim([0 xxlim]);
% ylim([-lim lim]);
% xlabel('$Dis.~\rm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\delta Inc~\rm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% title(name4);
% clear 'CE_record' ;'r2hill_record' ;strcat(tag,'_record_inout');
% 
