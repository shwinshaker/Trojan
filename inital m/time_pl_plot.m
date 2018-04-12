%%time_pl_plot
figure;
plot(Nep(:,1),mod(3*(Plu(:,5)+Plu(:,6)+Plu(:,7))-2*(Nep(:,5)+Nep(:,6)+Nep(:,7))-1*(Plu(:,5)+Plu(:,6)),360));
ylim([0 360]);set(gca,'ytick',0:10:360);
% figure;
% plot(Nep(:,1),mod(3*(elNep(:,5)+elNep(:,6)+elNep(:,7))-2*(Uran(:,5)+Uran(:,6)+Uran(:,7))-(elNep(:,5)+elNep(:,6)),360));
% ylim([0 360]);
% figure;
% plot(Nep(:,1),mod(3*(Nep(:,5)+Nep(:,6)+Nep(:,7))-2*(Jup(:,5)+Jup(:,6)+Jup(:,7))-(Nep(:,5)+Nep(:,6)),360));
% ylim([0 360]);