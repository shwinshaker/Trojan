clear;
ffname='RealPlutinosNpl';
% ffnameran='RanPlutinosNpl';
fname='1999CE119_1Gyr_40pl';
% fname='2001FU172_1Gyr_40pl';

di_file='di_record_inout.txt';
de_file='de_record_inout.txt';
de_file_ran='de_fit_perturb.txt';

%Dir='ServerMount';
Dir='swiftdata';

di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_file]);
%di_record_inout2=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname2,'/',di_file]);
de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',de_file]);
%de_record_inout_ran=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffnameran,'/',[fname,'_ran'],'/',de_file_ran]);
% de_record_inout_ran=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/','de_record_perturb.txt']);

fontsize=15;

%Gaussian

Nd=500;

%di_record_inout=abs(di_record_inout);
di_record_inout(abs(di_record_inout)>0.002)=[];
de_record_inout(abs(de_record_inout)>0.0001)=[];
% de_record_inout_ran(abs(de_record_inout_ran)>0.0001)=[];

%di_record_inout2(abs(di_record_inout2)>0.002)=[];

maxdi=max(abs(di_record_inout));
maxde=max(abs(de_record_inout));

dixx=0:Nd;
dixx=(-maxdi+dixx*2*maxdi/Nd);

dexx=0:Nd;
dexx=(-maxde+dexx*2*maxde/Nd);

%dixx=dixx*maxdi/Nd;

% xlim([-1 1]);
% ylim([0 yylim]);
countx=histcounts(di_record_inout,dixx)'/length(di_record_inout);
%countx2=histcounts(di_record_inout2,dixx)'/length(di_record_inout2);
countxe=histcounts(de_record_inout,dexx)'/length(de_record_inout);
% countxe_ran=histcounts(de_record_inout_ran,dexx)'/length(de_record_inout_ran);

dix=(dixx(2:end)+dixx(1:end-1))/2;
dex=(dexx(2:end)+dexx(1:end-1))/2;

figure;
set(gcf,'Position',[400,100,700,250],'color','w');


subplot(1,2,1);
plot(0,0,'w');
hold all;
bar(dix,countx,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
plot(dix,countx,'k-');
% ylim=get(gca,'ylim');
yylim=[0 0.03];
ylim(yylim);
plot([0 0],yylim,'k--');
xlabel('$\Delta I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');


subplot(1,2,2);
plot(0,0,'w');
hold all;
bar(dex,countxe,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
plot(dex,countxe,'k-');
% plot(dex,countxe_ran,'r-');
% ylim=get(gca,'ylim');
yylim=[0 0.03];
ylim(yylim);
plot([0 0],yylim,'k--');
xlabel('$\Delta e$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');


% plot(dix,countx2,'b.-');

% ylim([0 0.04]);
% %f = fittype('a*exp(-((x-b)/c)^2)','independent','x');
% f = fittype('a*x/(x^2)^2','independent','x');
% opt=fitoptions(f);
% opt.StartPoint=1;
% [cfun,gof] = fit(dix(:),countx(:),f,opt);
% %fitx = cfun.a*exp(-((dix-cfun.b)/cfun.c).^2);
% fitx = cfun.a*dix./(dix.^2).^2;
% 
% plot(dix,fitx,'r.');
hold off;
    
% for i=1:N
%     di_dis=di_record_perturb(dhill*(i-1) < dishill & dishill <= dhill*i);
%     %maxdi=max(abs(di_dis));
%     % the fit data of max di
%     %maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b+cfun_a.c)*maxmaxdi;
%     if strcmp(Boundary,'Infinity')
%         maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b)*maxmaxdi*MaxdiMult;
%     elseif strcmp(Boundary,'Threshold')
%         maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b-cfun_a.a./hillmax^cfun_a.b)*maxmaxdi*MaxdiMult;
%     elseif strcmp(Boundary,'None')
%         maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b+cfun_a.c)*maxmaxdi*MaxdiMult;
%     end
% 
%     disp(['maxdi:',num2str(maxdi)]);
%     if isempty(di_dis)
%         continue;
%     end
%     xx=1:Nd;
%     for j=1:Nd
%         countd(j)=length(find(-maxdi+(j-1)*2*maxdi/Nd < di_dis & di_dis <= -maxdi+j*2*maxdi/Nd));
%     end
%     x=(-maxdi+(xx-1/2)*2*maxdi/Nd)/maxdi;
%     xxitp=1:Nd*10;
%     xitp=(-maxdi+(xxitp-1/2)*2*maxdi/(Nd*10))/maxdi;
%     y=countd/length(di_dis);
%     [cfun,gof] = fit(x(:),y(:),f,opt);
%     yy = cfun.a*exp(-((xitp-cfun.b)/cfun.c).^2);
%     disp(cfun.a*cfun.c);
%     plot(x,y,'-','color',[(i-1)/N 0 1-(i-1)/N]);
% %     plot(xitp,yy,'-','color',[(i-1)/N 0 1-(i-1)/N]);
%     if cfun.b>=0
%         set(text(-0.9,(i+7)/(N+7)*yylim,['$$y=',num2str(cfun.a,'%.4f'),' \exp{\left\{-\left[(x-',num2str(cfun.b,'%.4f'),')/',num2str(abs(cfun.c),'%.4f'),'\right]^2\right\}}$$'])...,
%         ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N]);
%     else
%         set(text(-0.9,(i+7)/(N+7)*yylim,['$$y=',num2str(cfun.a,'%.4f'),' \exp{\left\{-\left[(x+',num2str(-cfun.b,'%.4f'),')/',num2str(abs(cfun.c),'%.4f'),'\right]^2\right\}}$$'])...,
%         ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N]);
%     end
%     set(text(0.9,(i+7)/(N+7)*yylim,['$$R^2=',num2str(gof.rsquare),'$$'])...,
%         ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N])
% end
% hold off;
