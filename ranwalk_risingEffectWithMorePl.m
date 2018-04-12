clear;

Np=power(10,0:0.005:6);

figure;
set(gcf,'Position',[400,100,700,350],'color','w');
fontsize=15;

mT=0.01;mPtot=1;
loglog(Np,(1+Np*mT/mPtot)/(1+mT/mPtot),'k-','linewidth',2);hold all;
loglog(Np,(1+Np*mT/mPtot).*Np.^(-1/2)/(1+mT/mPtot),'k--','linewidth',2);
xxlim=get(gca,'xlim');
plot(xxlim,[1 1],'k-.');

xlabel('$N_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$E(G_e)/E(G_e)_{N_P=1}$','fontsize',fontsize,'Interpreter','latex');
