%%bar eh
clear;
ffname='RealPlutinosNpl';
%ffname='RanPlutinos';
%ffname='fit_record_contrast_201606';
%ffname='fit_record_contrast';
%fname='2001KN77_1Gyr';
%fname='1999CE119_1Gyr';
%fname='1999CE119_1Gyr_5.0_hill';
%fname='1998HK151_1Gyr';
%fname='1999CE119_1Gyr_20pl';
fname='1999CE119_1Gyr_40pl';
%fname='1999CE119_40pl';
%fname='1999CE119_100k';
%fname='2001KN77_1Gyr_20pl';
CE_file='CE_record.txt';
di_file='di_record_perturb.txt';
%di_file='di_fit_perturb.txt';
%CE_file='ran_record.txt';
tag='real';%'ran';
%ran_init_file='ran_record.txt';

Dir='ServerMount';

CE_record=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',CE_file));
di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_file));

if strcmp(tag,'ran')
    ierror=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt'));
    CE_record(ierror,:)=[];
end

%ran_record=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',ran_init_file));
if strcmp(fname,'1999CE119_1Gyr_20pl')
    CE_record(54,:)=[];
end

% load(strcat('~/Documents/swiftdata/LAB/CE_realp/fit_record_contrast/',name,'/ran_CE_di.txt'));
r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);
hill=sqrt(r2hill_mean);
dishill=CE_record(:,2)/hill;
%dishill_init=ran_record(:,2)/hill; %%%% for ran specially

%% remove points too close
% dishillThresh=0.03;
% di_record_perturb(dishill<dishillThresh)=[];
% dishill(dishill<dishillThresh)=[];

%% di scale for the fit between distance and maxdi at this distance
maxmaxdi=0.05;
hillmax=3.5;
MaxdiMult=2.0; %% Scale for 3 fit = fitted maxdi in 2 * Multiplier
%Boundary='Infinity';
Boundary='Threshold';
fontsize=15;



figure(1);
set(gcf,'Position',[400,100,800,800],'color','w');

%yylim=0.025;
yylim=0.05;
N=50;
dhill=hillmax/N;
count=zeros(N,1);
%countraninit=zeros(N,1);
tot=length(dishill);
% countran=zeros(1,N);
for i=1:N
    count(i)=length(find(dhill*(i-1) < dishill & dishill <= dhill*i));
end
% for i=1:N
%     countraninit(i)=length(find(dhill*(i-1) < dishill_init & dishill_init <= dhill*i));
% end

% for i=1:N
%     countran(i)=length(find((i-1)/10*hill<ran_CE_di(:,1) & ran_CE_di(:,1)<=i/10*hill));
% end
x=1:N;
x=x';
y=count(x)/tot;
% yinit=countraninit(x)/tot;
x=x-1/2;
x=x*dhill;%/hillmax;
[P,H]=polyfit(x,y,1);
R=corrcoef(x,y);
yfit=polyval(P,x);
%semilogx(x/10,count(x),'k^');
% subplot(1,2,1);
% stem(x/10,count(x),'k.');

subplot(2,2,1);

% plot(x,yinit,'k.--');
%plot(x,(y+yinit)/2,'b.-');
bar(x,y,'facecolor',[0.9 0.9 0.9],'edgecolor','none');hold all;
plot(x,y,'k.--','linewidth',1.3);
plot(x,yfit,'r-');
hold off;
%axis square;
%xlim([0.1 3.5]);
ylim([0 yylim]);
%axis([0 N 0 max(count)]);
%ylim=get(gca,'ylim');
%yylim=ylim(2);
% set(gca,'xticklabel',0:estep:N*estep);
set(text(0.3,8/9*yylim,['$$N_{CE}= ',num2str(length(dishill)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
if P(2)>=0
    set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x+',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
else
    set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
end
set(text(0.3,3.5/5*yylim,['$$R^2=',num2str(R(1,2)*R(2,1),'%.4f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%title(num2str(length(CErecord),'N=%04d'));
xlabel('$$\mathrm{Dis.~(R_H)}$$','fontsize',fontsize,'Interpreter','latex');
ylabel('Proportion','fontsize',fontsize,'Interpreter','latex');


%% max di fit
subplot(2,2,2);

% power fit !! Force the value at hillmax equals to 0
if strcmp(Boundary,'Infinity')
    fstr='a/x^b'; %% Boundary condition: Infinity = 0   
elseif strcmp(Boundary,'Threshold')
    fstr=['a/x^b-a/',num2str(hillmax),'^b'];
elseif strcmp(Boundary,'None')
    fstr='a/x^b+c';
end

f = fittype(fstr,'independent','x');
opt=fitoptions(f);
opt.StartPoint=[1 1];
N=30; %% distance sample points
yylim=1;

dhill=hillmax/N;
%maxmaxdi=max(abs(di_record_perturb));
maxdi_posi=zeros(N,1);
maxdi_nega=zeros(N,1);
maxdi_abs=zeros(N,1);

%% pre plot
plot(0,0,'w');axis square;
hold all;
xlim([0 hillmax]);
ylim([-yylim yylim]);
plot([0 hillmax],[0 0],'k--');
xlabel('Dis /Hill','Interpreter','latex','fontsize',fontsize);
ylabel(['$$\rm{|di|_{max}/',num2str(maxmaxdi),'DEG}$$'],'Interpreter','latex','fontsize',fontsize);

for i=1:N
    di_dis=di_record_perturb(dhill*(i-1) < dishill & dishill <= dhill*i);
    di_posi=di_dis(di_dis>=0);
    di_nega=di_dis(di_dis<0);
    if isempty(di_posi)
        maxdi_posi(i)=0.0;
    else
        maxdi_posi(i)=max(di_posi);
    end
    if isempty(di_nega)
        maxdi_nega(i)=0.0;
    else
        maxdi_nega(i)=-max(-di_nega);
    end
    if isempty(di_dis)
        maxdi_abs(i)=0.0;
    else
        maxdi_abs(i)=max(abs(di_dis));
    end
end
xx=1:N;
xx=xx-1/2;
x=dhill*xx;
y_posi=maxdi_posi/maxmaxdi;
y_nega=-maxdi_nega/maxmaxdi;
y_abs=maxdi_abs/maxmaxdi;
[cfun_p,gof_p] = fit(x(:),y_posi(:),f,opt);
[cfun_n,gof_n] = fit(x(:),y_nega(:),f,opt);
[cfun_a,gof_a] = fit(x(:),y_abs(:),f,opt);
xxitp=1:N*10;
xxitp=xxitp-1/2;
xitp=hillmax/(N*10)*xxitp;
if strcmp(Boundary,'Infinity')
    yy_p = cfun_p.a./xitp.^cfun_p.b;
    yy_n = cfun_n.a./xitp.^cfun_n.b;
    yy_a = cfun_a.a./xitp.^cfun_a.b; 
elseif strcmp(Boundary,'Threshold')
    yy_p = cfun_p.a./xitp.^cfun_p.b-cfun_p.a./hillmax.^cfun_p.b;
    yy_n = cfun_n.a./xitp.^cfun_n.b-cfun_n.a./hillmax.^cfun_n.b;
    yy_a = cfun_a.a./xitp.^cfun_a.b-cfun_a.a./hillmax.^cfun_a.b; 
elseif strcmp(Boundary,'None')
    yy_p = cfun_p.a./xitp.^cfun_p.b+cfun_p.c;
    yy_n = cfun_n.a./xitp.^cfun_n.b+cfun_n.c;
    yy_a = cfun_a.a./xitp.^cfun_a.b+cfun_a.c; 
end

plot(x,y_posi,'k.');
plot(x,-y_nega,'k.');
plot(x,y_abs,'k^');
plot(xitp,yy_p,'r-');
plot(xitp,-yy_n,'r-');
plot(xitp,yy_a,'r--');


if strcmp(Boundary,'Infinity')
    set(text(0.4,7/8*yylim,['$$y = ',num2str(cfun_p.a),' / x^{',num2str(cfun_p.b),'}$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,-7/8*yylim,['$$y = -[',num2str(cfun_n.a),' / x^{',num2str(cfun_n.b),'}]$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,5/8*yylim,['$$y = ',num2str(cfun_a.a),' / x^{',num2str(cfun_a.b),'}$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
elseif strcmp(Boundary,'Threshold')
    set(text(0.4,7/8*yylim,['$$y = ',num2str(cfun_p.a),' / x^{',num2str(cfun_p.b),'}+Const$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,-7/8*yylim,['$$y = -[',num2str(cfun_n.a),' / x^{',num2str(cfun_n.b),'}+Const]$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,5/8*yylim,['$$y = ',num2str(cfun_a.a),' / x^{',num2str(cfun_a.b),'}+Const$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
elseif strcmp(Boundary,'None')
    set(text(0.4,7/8*yylim,['$$y = ',num2str(cfun_p.a),' / x^{',num2str(cfun_p.b),'}',num2str(cfun_p.c),'$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,-7/8*yylim,['$$y = -[',num2str(cfun_n.a),' / x^{',num2str(cfun_n.b),'}',num2str(cfun_n.c),']$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
    set(text(0.4,5/8*yylim,['$$y = ',num2str(cfun_a.a),' / x^{',num2str(cfun_a.b),'}',num2str(cfun_a.c),'$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
end




set(text(0.4,6/8*yylim,['$$R^2=',num2str(gof_p.rsquare),'$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
set(text(0.4,-6/8*yylim,['$$R^2=',num2str(gof_n.rsquare),'$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
set(text(0.4,4/8*yylim,['$$R^2=',num2str(gof_a.rsquare),'$$']),'Interpreter','latex',...,
        'fontsize',fontsize,'color','r');
hold off;


%% di bar at some distance
subplot(2,2,3);
plot(0,0,'w');axis square;
hold all;
%Gaussian
f = fittype('a*exp(-((x-b)/c)^2)','independent','x');
opt=fitoptions(f);
opt.StartPoint=[1 0 1];

N=5; %% distance sample points
Nd=50;

countd=zeros(Nd,1);
yylim=0.5;

dhill=hillmax/N;

xlim([-1 1]);
ylim([0 yylim]);
xlabel('di/max|di|','fontsize',fontsize);
ylabel('CE times / Tot times within this dis. range','fontsize',fontsize);

for i=1:N
    di_dis=di_record_perturb(dhill*(i-1) < dishill & dishill <= dhill*i);
    %maxdi=max(abs(di_dis));
    % the fit data of max di
    %maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b+cfun_a.c)*maxmaxdi;
    if strcmp(Boundary,'Infinity')
        maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b)*maxmaxdi*MaxdiMult;
    elseif strcmp(Boundary,'Threshold')
        maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b-cfun_a.a./hillmax^cfun_a.b)*maxmaxdi*MaxdiMult;
    elseif strcmp(Boundary,'None')
        maxdi = (cfun_a.a./(dhill*(i-1/2)).^cfun_a.b+cfun_a.c)*maxmaxdi*MaxdiMult;
    end

    disp(['maxdi:',num2str(maxdi)]);
    if isempty(di_dis)
        continue;
    end
    xx=1:Nd;
    for j=1:Nd
        countd(j)=length(find(-maxdi+(j-1)*2*maxdi/Nd < di_dis & di_dis <= -maxdi+j*2*maxdi/Nd));
    end
    x=(-maxdi+(xx-1/2)*2*maxdi/Nd)/maxdi;
    xxitp=1:Nd*10;
    xitp=(-maxdi+(xxitp-1/2)*2*maxdi/(Nd*10))/maxdi;
    y=countd/length(di_dis);
    [cfun,gof] = fit(x(:),y(:),f,opt);
    yy = cfun.a*exp(-((xitp-cfun.b)/cfun.c).^2);
    disp(cfun.a*cfun.c);
    plot(x,y,'-','color',[(i-1)/N 0 1-(i-1)/N]);
%     plot(xitp,yy,'-','color',[(i-1)/N 0 1-(i-1)/N]);
    if cfun.b>=0
        set(text(-0.9,(i+7)/(N+7)*yylim,['$$y=',num2str(cfun.a,'%.4f'),' \exp{\left\{-\left[(x-',num2str(cfun.b,'%.4f'),')/',num2str(abs(cfun.c),'%.4f'),'\right]^2\right\}}$$'])...,
        ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N]);
    else
        set(text(-0.9,(i+7)/(N+7)*yylim,['$$y=',num2str(cfun.a,'%.4f'),' \exp{\left\{-\left[(x+',num2str(-cfun.b,'%.4f'),')/',num2str(abs(cfun.c),'%.4f'),'\right]^2\right\}}$$'])...,
        ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N]);
    end
    set(text(0.9,(i+7)/(N+7)*yylim,['$$R^2=',num2str(gof.rsquare),'$$'])...,
        ,'Interpreter','latex','fontsize',fontsize-4,'color',[(i-1)/N 0 1-(i-1)/N])
end
hold off;

% subplot(2,2,4);
% plot(0,0,'w');axis square;
% hold all;
% maxdi=max(abs(di_record_perturb));
% N=5;
% Nbase=N/3.5;
% Nd=100;
% countd=zeros(Nd,1);
% minlog=log(min(abs(di_record_perturb))/maxdi);
% for i=1:N
%     di_dis=di_record_perturb(dhill*(i-1) < dishill & dishill <= dhill*i);
%     if isempty(di_dis)
%         continue;
%     end
%     xx=1:Nd;
%     for j=1:Nd
%         countd(j)=length(find(minlog-(j-1)*minlog/Nd < log(abs(di_dis)/maxdi) & log(abs(di_dis)/maxdi) <= minlog-j*minlog/Nd));
%     end
%     plot(minlog-(xx-1/2)*minlog/Nd,countd/length(di_dis),'-','color',[(i-1)/N 0 1-(i-1)/N]);
% end
% hold off;
% xlabel('log( |di| /max |di| )','fontsize',fontsize);
% ylabel('CE times / Tot times within this dis. range','fontsize',fontsize);

% subplot(1,2,2);
% stem(x/10,countran(x),'k.');
% axis square;
% xlim([0.1 3.5]);
% set(gca,'ylim',ylim_1);
% set(text(0.5,8/9*(max(count(x))),strcat('N_{CE}= ',num2str(length(ran_CE_di)))),'color','red');
% title('ran-CE-bar');
