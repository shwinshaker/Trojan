%% plot MI

x=exp(log(1e-3):0.1:log(1e3));

fontsize=15;
figure;
set(gcf,'Position',[400,100,700,500],'color','w');
h1=loglog(x,x.^2.*(6.8+log(x)/3*2),'k-');hold all;
h2=loglog(x,x.^2*6.8,'k--');
h3=loglog(x,x.^2.1*6.8,'r-');
hold off;
xlabel('$\nu_P$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_I/\chi_M$','fontsize',fontsize,'Interpreter','latex');
yylim=get(gca,'ylim');
set(text(xp,yylim(2)/10,['$$(1)~~\frac{M_I}{\chi_M}={\nu_P}^2\left(6.8+\frac{2}{3}\ln{\nu_P}\right)$$']),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,yylim(2)/10^3,'$$(2)~~\frac{M_I}{\chi_M}=6.8{\nu_P}^2$$'),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');
set(text(xp,yylim(2)/10^5,'$$(3)~~\frac{M_I}{\chi_M}=6.8{\nu_P}^{2.1}$$'),...,
    'Interpreter','latex','fontsize',fontsize,'color','red');

legend([h1 h2 h3],{'(1)','(2)','(3)'},'fontsize',fontsize,'Interpreter','latex','location','southeast')