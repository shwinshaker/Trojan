
clear;

ffname='symba_intro_fictiousTrojans_fast';
fname={'NoPlutino_1Gyr','Pluto_10Mass_1Gyr','Pluto_100Mass_1Gyr', ...,
    'Multi-1+9Plutino_100M_1Gyr', 'Multi-1+9Plutino_10M_1Gyr'};
for i=1:length(fname)
    meanEle=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{i},'/mean.txt']);
    eval(['meanEle',int2str(i),'=meanEle;']);
end

figure;
set(gcf,'Position',[400,100,700,500],'color','w');
fontsize=15;
colorspec={'k','b','r','m','c'};

subplot(2,2,1)
plot(mean(meanEle(:,4)), mean(meanEle(:,5)),'w+');hold all;
for i=1:length(fname)
    eval(['meanEle=meanEle',int2str(i),';']);
    h=plot(meanEle(:,4),meanEle(:,5),[colorspec{i},'.']);
    eval(['h',int2str(i),'=h;']);
end
hold off;
xlabel('$a~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$e$','fontsize',fontsize,'Interpreter','latex');
legend([h1 h2 h3 h4 h5],{'$S+N$','$S+N+P(10)$','$S+N+P(100)$','$S+N+10P(10)$','$S+N+10P(1)$'}, ...,
        'fontsize',12,'location','southeast','box','on','Interpreter','latex');


subplot(2,2,2)
plot(mean(meanEle(:,4)), mean(meanEle(:,6)),'w+');hold all;
for i=1:length(fname)
    eval(['meanEle=meanEle',int2str(i),';']);
    h=plot(meanEle(:,4),meanEle(:,6),[colorspec{i},'.']);
    eval(['h',int2str(i),'=h;']);
end
hold off;
xlabel('$a~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$inc~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');


subplot(2,2,3)
% plot(mean(meanEle(:,5)), mean(meanEle(:,6)),'w+');hold all;
meane=[meanEle1(:,5);meanEle2(:,5);meanEle3(:,5);meanEle4(:,5);meanEle5(:,5)];
meani=[meanEle1(:,6);meanEle2(:,6);meanEle3(:,6);meanEle4(:,6);meanEle5(:,6)];
groupind=cell(length(meanEle1)+length(meanEle2)+length(meanEle3)+length(meanEle4)+length(meanEle5),1);
tag={'S+N','S+N+P(10)','S+N+P(100)','S+N+10P(10)','S+N+10P(1)'};
j0=1;
for i=1:length(fname)
    eval(['meanEle=meanEle',int2str(i),';']);
    for j=j0:j0+length(meanEle)-1
        groupind{j}=tag{i};
    end
    j0=j+1;
end
% groupind=[ones(length(meanEle1),1); ones(length(meanEle2),1)*2; ones(length(meanEle3),1)*3; ...,
%     ones(length(meanEle4),1)*4; ones(length(meanEle5),1)*5];

scatterhist(meane,meani,'Group',groupind,'Location','SouthWest',...
    'Direction','out','Color',colorspec,'LineStyle',{'-','-','-.','-.','-'},...
    'LineWidth',ones(5,1),'Marker','.+xx+','MarkerSize',ones(5,1)*7);
% for i=1:length(fname)
%     eval(['meanEle=meanEle',int2str(i),';']);
%     h=plot(meanEle(:,5),meanEle(:,6),[colorspec{i},'.']);
%     eval(['h',int2str(i),'=h;']);
% end
% hold off;
xlabel('$e$','fontsize',fontsize,'Interpreter','latex');
ylabel('$inc~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');


% % sanity check
% if ~(all(mean1(:,1)==mean2(:,1)) && all(mean1(:,1)==mean3(:,1)))
%     error('id not consistent!')
% end
% % 
% idCmpa=zeros(length(meanEle),4);
% idCmpa(:,1)=meanEle1(:,1);
% idCmpa(:,2)=meanEle2(:,4);
% idCmpa(:,3)=meanEle3(:,4);
% idCmpa(:,4)=meanEle4(:,4);
% 
% idCmpe=zeros(length(meanEle),4);
% idCmpe(:,1)=meanEle1(:,1);
% idCmpe(:,2)=meanEle2(:,5);
% idCmpe(:,3)=meanEle3(:,5);
% idCmpe(:,4)=meanEle4(:,5);
% 
% idCmpi=zeros(length(meanEle),4);
% idCmpi(:,1)=meanEle1(:,1);
% idCmpi(:,2)=meanEle1(:,6);
% idCmpi(:,3)=meanEle2(:,6);
% idCmpi(:,4)=meanEle3(:,6);
% 
% % idCmpa=sortrows(idCmpa,2);
% % idCmpe=sortrows(idCmpe,2);
% % idCmpi=sortrows(idCmpi,2);
% % 
% figure;
% set(gcf,'Position',[400,100,700,500],'color','w');
% 
% subplot(2,2,1)
% plot(idCmpa(:,2),'b.');hold on;
% plot(idCmpa(:,3),'g.');
% plot(idCmpa(:,4),'r.');
% 
% subplot(2,2,2)
% plot(idCmpe(:,2),'b.');hold on;
% plot(idCmpe(:,3),'g.');
% plot(idCmpe(:,4),'r.');
% 
% subplot(2,2,3)
% plot(idCmpi(:,2),'b.');hold on;
% plot(idCmpi(:,3),'g.');
% plot(idCmpi(:,4),'r.');
% 
