% %time_plot
% 
% figure;
% subplot(4,1,1);
% fl=2;
% plot(el(:,1),el(:,fl),'r',elNep(:,1),elNep(:,fl),'k')...,
% %     ,elr(:,1),elr(:,fl), ...,
% % 'b',elno(:,1),elno(:,fl),'g',elNepno(:,1),elNepno(:,fl),'k.')
% ylim([29 31]);
% subplot(4,1,2);
% fl=3;
% plot(el(:,1),el(:,fl),'r',elNep(:,1),elNep(:,fl),'k') ...,
% %     ,elr(:,1),elr(:,fl), ...,
% % 'b',elno(:,1),elno(:,fl),'g',elNepno(:,1),elNepno(:,fl),'k.')
% ylim([0 0.1]);
% subplot(4,1,3);
% fl=4;
% plot(el(:,1),el(:,fl),'r',elNep(:,1),elNep(:,fl),'k') ...,
% %     ,elr(:,1),elr(:,fl), ...,
% % 'b',elno(:,1),elno(:,fl),'g',elNepno(:,1),elNepno(:,fl),'k.')
% ylim([0 20]);hold on;
% plot(get(gca,'xlim'),[4 4]);
% subplot(4,1,4)
% plot(el(:,1),mod(el(:,5)+el(:,6)+el(:,7) ...,
%     -elNep(:,5)-elNep(:,6)-elNep(:,7),360),'r');hold on;ylim([0 360]);
% % plot(elr(:,1),mod(elr(:,5)+elr(:,6)+elr(:,7) ...,
% %     -elNep(:,5)-elNep(:,6)-elNep(:,7),360),'b');hold on;
% % plot(elno(:,1),mod(elno(:,5)+elno(:,6)+elno(:,7) ...,
% %     -elNepno(:,5)-elNepno(:,6)-elNepno(:,7),360),'k');hold on;
% plot(get(gca,'xlim'),[300 300],'k');hold on;
% plot(get(gca,'xlim'),[60 60],'k');hold on;
% % plot(get(gca,'xlim'),[-360 -360],'k');

%%plot Nep and random plutinos
subplot(6,1,1);
plot(Nep(:,1),Nep(:,2));ylim([0 50]);
subplot(6,1,2);
plot(plutino1(:,1),plutino1(:,2));ylim([0 50]);
subplot(6,1,3);
plot(plutino2(:,1),plutino2(:,2));ylim([0 50]);
subplot(6,1,4);
plot(plutino3(:,1),plutino3(:,2));ylim([0 50]);
subplot(6,1,5);
plot(plutino4(:,1),plutino4(:,2));ylim([0 50]);
subplot(6,1,6);
plot(plutino5(:,1),plutino5(:,2));ylim([0 50]);


% syms nplutino;
% nplutino=5;
% for i=2,nplutino;
% eval(['subplot(6,1,',num2str(i),')']);
% eval(['plot(plutino',num2str(i-1),'(:,1),plutino',num2str(i-1),'(:,2))']);
% ylim([0 50]);hold on;
% end