%%temp
plot(mA4(:,1),mA4(:,2),'k');ylim([0 5e4]);xlabel('t/s');ylabel('ÊýÖµ');
Rp=(sum(mA4(:,2))-sum(mA0(:,2)))/50;
Pi=3.96e-19*Rp/0.15;
Nt=sum(mA4(:,2));
Nd=sum(mA0(:,2));
snr=(Nt-Nd)/sqrt(Nt+Nd);