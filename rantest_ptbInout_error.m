clear;

ffname='RealPlutinosNpl';
fname='1999CE119_1Gyr_40pl';
%fname='2001FU172_1Gyr_40pl';

tag='di';

CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);
hill=sqrt(r2hill_mean);

ierror=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt']);
CE_record(ierror,:)=[];

di_record_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt']);
di_fit_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt']);
di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt']);

di_record_perturb(ierror)=[];
di_record_inout(ierror)=[];

tag='de';

de_record_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt']);
de_fit_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt']);
de_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt']);

de_record_perturb(ierror)=[];
de_record_inout(ierror)=[];

fontsize=15;
xxlim=3.5;
yylim=[1e-5 1e3];
Row=1;
Column=2;

figure(1);
set(gcf,'Position',[400,100,700,250],'color','w');
% annotation(gcf,'textbox','String',{[fname,' ',tag]},'FontSize',12,'Position',[0.35 0.88 0.10 0.10],'edgecolor',get(gcf,'color'))

subplot(Row,Column,1);

rlterrdi=abs((di_fit_perturb-di_record_inout)./di_record_inout);

semilogy(CE_record(:,2)/hill,rlterrdi,'k.');hold on;
rlterrsortdi=sort(rlterrdi);
markdi=rlterrsortdi(floor(length(rlterrsortdi)*0.90));
plot([0 xxlim],[markdi markdi],'r-','linewidth',3);

xlim([0 xxlim]);
ylim(yylim);
xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\eta_I$','fontsize',fontsize,'Interpreter','latex');

subplot(Row,Column,2);

rlterrde=abs((de_fit_perturb-de_record_inout)./de_record_inout);

semilogy(CE_record(:,2)/hill,rlterrde,'k.');hold on;
rlterrsortde=sort(rlterrde);
markde=rlterrsortde(floor(length(rlterrsortde)*0.90));
plot([0 xxlim],[markde markde],'r-','linewidth',3);

xlim([0 xxlim]);
ylim(yylim);
xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\eta_e$','fontsize',fontsize,'Interpreter','latex');

figure;
set(gcf,'Position',[400,100,700,250],'color','w');
subplot(Row,Column,1);
plot(CE_record(:,2)/hill,(di_fit_perturb),'k.');hold on;
subplot(Row,Column,2);

plot(CE_record(:,2)/hill,(di_record_inout),'r.');

% subplot(Row,Column,2);
% semilogy(CE_record(:,2)/hill,abs(de_fit_perturb),'k.');hold on;
% semilogy(CE_record(:,2)/hill,abs(de_record_inout),'r.');

%loglog(di_record_inout,di_record_perturb./di_record_inout,'k.');
% 
% pick=de_record_perturb~=de_record_inout;
% de_record_perturb=de_record_perturb(pick);
% de_record_inout=de_record_inout(pick);
% 
% pick0=de_record_inout~=0;
% de_record_inout=de_record_inout(pick0);
% de_record_perturb=de_record_perturb(pick0);
% 
% errde=((de_record_perturb-de_record_inout)./de_record_inout);
% errdi=((di_record_perturb-di_record_inout)./di_record_inout);
% 
% eta=median(errde)/median(errdi);
% errde_norm=errde/eta;
% 
% de_record_perturb_revise=errde_norm.*de_record_inout+de_record_inout;
% 
% plot(de_record_inout,de_record_perturb_revise,'k.');hold on;
% plot(de_record_inout,de_record_perturb,'r.');

% err_revise=(de_record_perturb_revise-de_record_perturb)./de_record_perturb;

% pick0=de_record_perturb>0;
% de_record_inout=de_record_inout(pick0);
% de_record_perturb=de_record_perturb(pick0);
% 
% pick=de_record_perturb~=de_record_inout;
% de_record_perturb=de_record_perturb(pick);
% de_record_inout=de_record_inout(pick);
% 
% err_init=(de_record_perturb-de_record_inout)./de_record_inout;
% 
% err_revise=(de_record_inout-de_record_perturb)./de_record_perturb;
% 
% [P,H]=polyfit(log(abs(de_record_perturb)),log(abs(err_revise)),1);  
% %R=corrcoef(log(Fitx),log(TimesListFit));
% err_revise_fit=exp(polyval(P,log(abs(de_record_perturb))));
% 
% de_record_perturb_revise=de_record_perturb.*(1+err_revise_fit*rand());
% % plot(de_record_inout,de_record_perturb_revise,'k.');hold on;
% % plot(de_record_inout,de_record_perturb,'r.');
% err_fit=(de_record_perturb_revise-de_record_inout)./de_record_inout;
% loglog(de_record_inout,err_fit,'k.');hold on;
% loglog(de_record_inout,err_init,'r.');
%plot(de_record_perturb,err_revise,'k.');
% log(abs(de_record_perturb),abs(err_revise),'k.');hold on;
% log(abs(de_record_perturb),abs(err_revise_fit),'r.');

%loglog(abs(di_record_inout)/max(abs(di_record_inout)),errdi,'c.');

% loglog(abs(de_record_inout)/max(abs(de_record_inout)),errde,'k.');hold on;
% loglog(abs(di_record_inout)/max(abs(di_record_inout)),errdi,'c.');

% loglog(abs(de_record_inout)/max(abs(de_record_inout)),abs(de_record_perturb./de_record_inout),'k.');hold on;
% loglog(abs(di_record_inout)/max(abs(di_record_inout)),abs(di_record_perturb./di_record_inout),'c.');

% semilogy(CE_record(:,2)/hill,di_record_perturb./di_record_inout,'k.');
% % set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval([tag,'_fit_inout']))))),'color','red');
% % set(text(1.5,5/6*lim,strcat('N=',num2str(length(CE_record)))),'color','red');
% hold on;
% plot([0 xxlim],[1 1],'r-','linewidth',2);
% xlim([0 xxlim]);
% ylim(yylim);
% xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\eta$','fontsize',fontsize,'Interpreter','latex');

%title('$\Delta i$','fontsize',fontsize,'Interpreter','latex');

% subplot(Row,Column,2);
% semilogy(CE_record(:,2)/hill,abs((de_record_perturb-de_record_inout)./de_record_inout),'k.');
% % set(text(1.5,3/4*lim,strcat('sum=',num2str(sum(eval([tag,'_fit_inout']))))),'color','red');
% % set(text(1.5,5/6*lim,strcat('N=',num2str(length(CE_record)))),'color','red');
% 
% %grid on;
% xlim([0 xxlim]);
% ylim(yylim);
% xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$Err$','fontsize',fontsize,'Interpreter','latex');
% title('$\Delta e$','fontsize',fontsize,'Interpreter','latex');


% maxabsdi=max(abs(de_record_inout));
% 
% N=500;
% di_norm_0=de_record_inout(de_record_inout~=0);
% di_norm_02=de_record_perturb(de_record_perturb~=0);
% Minlog=-log(min(min(abs(di_norm_0)),min(abs(di_norm_02))));
% 
% dix=0:N;
% dix=dix';
% dlog=Minlog/N;
% dix=exp(-dix*dlog);
% dix=sort(dix);
% dixx=(dix(2:end)+dix(1:end-1))/2;
% % cdf(de_record_inout,dixx);
% 
% countiot=histcounts(abs(di_norm_0),dix)'/length(de_record_inout);
% countptb=histcounts(abs(di_norm_02),dix)'/length(de_record_perturb);
% semilogx(dixx,countiot,'k.-','linewidth',1.3,'markersize',10);hold on;
% semilogx(dixx,countptb,'r.-','linewidth',1.3,'markersize',10);
