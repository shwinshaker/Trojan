%% total i change
clear;

ffname='RealPlutinosNpl';
fname='2001FU172_1Gyr_40pl';
tag='di';

CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
ddata=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/','di_record_perturb.txt'));

if strcmp(tag,'di')
%     limdi=0.5;%0.03;
%     limsum=1.0;%0.2;
%     limDi=0.5;
    limdi=0.1;
    limsum=0.2;
    limDi=0.1;
    label='I';
    unit=' (DEG)';
else
    limdi=0.001;
    limDi=0.001;
    limsum=0.02;
    label='e';
    unit='';
end

fontsize=15;
N=200;

time=CE_record(:,1)/365.25; %% time scale
maxtime=1e9;
di_sum=cumsum(ddata);

%% inc total change in each time interval

Di=zeros(1,N);
for i=1:N
    Di(i)=sum(ddata((i-1)/N*maxtime<time & time<=i/N*maxtime));
end

BottomRetainWidth=0.05;
LeftRetainWidth=0.2;
Height=0.23;
Width=0.7;

figure(2);
set(gcf,'Position',[400,100,700,400],'color','w');

row=3;

subplot(row,1,1);

plot(time/maxtime,di_sum,'k-');hold on;
plot([0 1],[0 0],'k--');

% ylim([-limsum limsum]);
%ylabel(['$A_',label,'=\sum_{i=1}^k \Delta ',label,'_i~\mathrm{',unit,'}$'],'fontsize',fontsize,'Interpreter','latex');
ylabel(['$A_{',label,',n}','~\mathrm{',unit,'}$'],'fontsize',fontsize,'Interpreter','latex');

yylim=get(gca,'ylim');
[maxsum,imax]=max(di_sum);
[minsum,imin]=min(di_sum);
plot([time(imax)/maxtime-0.05 time(imax)/maxtime+0.05],[maxsum maxsum],'r-');
plot([time(imax)/maxtime-0.05 time(imin)/maxtime+0.05],[minsum minsum],'r-');
xposition=[LeftRetainWidth+(time(imax)/maxtime)*Width LeftRetainWidth+(time(imax)/maxtime)*Width];
yposition=[(maxsum-yylim(1))/(yylim(2)-yylim(1))*Height+BottomRetainWidth+3*Height (minsum-yylim(1))/(yylim(2)-yylim(1))*Height+BottomRetainWidth+3*Height];

annotation('doublearrow',xposition,yposition,'color','r');
%text(imax/length(di_sum),maxsum,' \arrow sin(\pi)');
text(time(imax)/maxtime+0.005,0,'$G_I$','Interpreter','latex','color','r');

set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth*0.8 0.03+3*Height 0.05 0.05],'edgecolor','none','string',...,
    '(A)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

ytick=get(gca,'ytick');
set(gca,'yTick',ytick(2:end-1));


subplot(row,1,2); %% bar is better?
Ni=1:N;
Ni=Ni-1/2;
% plot(Ni/N,Di,'b-');hold on;
bar(Ni/N,Di,'k');hold on;
plot([0 1],[0 0],'k-');

% ylim([-limDi limDi]);
ylabel(['$\Xi ',label,'~\mathrm{',unit,'}$'],'fontsize',fontsize,'Interpreter','latex','Interpreter','latex');

set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth*0.8 0.03+2*Height 0.05 0.05],'edgecolor','none','string',...,
   '(B)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

ytick=get(gca,'ytick');
set(gca,'yTick',ytick(2:end-1));

subplot(row,1,3);
plot(time/maxtime,ddata,'k.');hold on;
plot([0 1],[0 0],'k-');
% ylim([-limdi limdi]);
xlabel('$t/t_{max}$','fontsize',fontsize,'Interpreter','latex');
ylabel(['$\Delta ',label,'~\mathrm{',unit,'}$'],'fontsize',fontsize,'Interpreter','latex','Interpreter','latex');

set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
annotation('textbox',[Width+LeftRetainWidth*0.8 0.03+Height 0.05 0.05],'edgecolor','none','string',...,
    '(C)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

ytick=get(gca,'ytick');
set(gca,'yTick',ytick(2:end-1));

