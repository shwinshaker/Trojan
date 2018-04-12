figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
fontsize=15;
subplot(2,1,1);
plot(ePL(:,2),IPL(:,2)/pi*180,'k.');hold on;plot(EccList,30.943,'k+')
yylim=[20 32];
ylim(yylim);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
set(gca,'xticklabel',[]);
set(gca,'yTick',yylim(1)+2:2:yylim(2));
ylabel('$i~\rm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');

subplot(2,1,2);
plot(ePL(:,2),aPL(:,2),'k.');hold on;plot(EccList,39.636,'k+')
yylim=[39.4 39.7];
ylim(yylim);
set(gca,'yTick',yylim(1):0.05:yylim(2)-0.05);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
ylabel('$a~\rm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
xlabel('$e$','fontsize',fontsize,'Interpreter','latex');