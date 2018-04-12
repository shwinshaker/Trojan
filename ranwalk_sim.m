
if ~(exist('di_record_inout','var') && exist('CE_record','var'))
    
    clear;

    %ffname='SimRanWalk';
    ffname='RanSimRanWalk';
    fname='1999CE119_2004UP10_da_2';
    %dirname='ServerMount/LAB/CE_realp';
    dirname='swiftdata/Trojan/LAB/CE_realp';
%     ffname2='RanPlutinosNK';
%     %fname2='1999CE119_2004UP10_3MP';
%     fname2='1999CE119_1M';
%     fname3='1999CE119_1Gyr_40pl';
    diname='di_record_inout';
    dename='de_record_inout';
    di_record_inout=load(['~/Documents/',dirname,'/',ffname,'/',fname,'/',diname,'.txt']);
    CE_record=load(['~/Documents/',dirname,'/',ffname,'/',fname,'/CE_record.txt']);
    tpel=load(['~/Documents/',dirname,'/',ffname,'/',fname,'/tpel.txt']);
    Nepel=load(['~/Documents/',dirname,'/',ffname,'/',fname,'/Nepel.txt']);
    
%     di_record_inout2=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname2,'/di_fit_perturb.txt']);
%     de_record_inout2=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname2,'/de_fit_perturb.txt']);
%     
%     de_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/RealPlutinosNpl/',fname3,'/',dename,'.txt']);

end

%% handle data
len=min(length(di_record_inout),length(CE_record));

di=di_record_inout(1:len,1);
de=di_record_inout(1:len,2);
da=di_record_inout(1:len,3);

%de2=di_record_inout2(:,3);

%time=di_record_inout(:,1)/3.6525e2;
time=CE_record(1:len,1)/3.6525e2;

a=tpel(:,2)-30;
e=tpel(:,3);
inc=tpel(:,4)/180*pi;
etime=tpel(:,1)/3.6525e2;

phi=mod((tpel(:,5)+tpel(:,6)+tpel(:,7))-(Nepel(:,5)+Nepel(:,6)+Nepel(:,7)),360);

% %% Plot de distribution
% figure;
% % plot(time,cumsum(di),'k.');hold all;
% % plot(time,cumsum(de),'r.');
% % plotyy(time,cumsum(di),time,cumsum(de));
% 
% N=100;
% 
% absde=abs(de);
% Minlog=-log(min(abs(absde)));
% 
% dix=0:N;
% dix=dix';
% dlog=Minlog/N;
% dix=exp(-dix*dlog);
% dix=sort(dix);
% dixx=(dix(2:end)+dix(1:end-1))/2;
% 
% countx=histcounts(absde,dix)'/length(absde);
% semilogx(dixx,countx,'k.-');
% 
% hold all;
% 
% %% RealPlutinosNplContrast
% absde2=abs(de_record_inout);
% countx2=histcounts(absde2,dix)'/length(absde2);
% semilogx(dixx,countx2,'r.-');
% 
% absde3=abs(de_record_inout2);
% countx3=histcounts(absde3,dix)'/length(absde3);
% semilogx(dixx,countx3,'k.-');

%% Contrast of positive and negative branch of sim CE
% dep=de_record_inout2(de_record_inout2>0);
% den=-de_record_inout2(de_record_inout2<0);
% countxp=histcounts(dep,dix)'/length(dep);
% countxn=histcounts(den,dix)'/length(den);
% semilogx(dixx,countxp,'r.-');hold on;
% semilogx(dixx,countxn,'b.-');
% 

%% fetch the CE times

temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
if isempty(temp)
    ejectNo=size(tpel,1);
else
    ejectNo=temp;
end
clear temp;
ejecttime=etime(ejectNo);
NCE=find(time<ejecttime,1,'last');

    
%% Plot random walk and el change
figure;
set(gcf,'Position',[400,100,700,800],'color','w');
fontsize=15;

subplot(3,1,1);
plot(time,cumsum(da),'k-');hold all;
plot(etime,a,'g-');
% Max N da values
[val,ind]=sort(abs(da),'descend');
Nmax=10;
ploty=zeros(1,Nmax*3);
plotx=zeros(1,Nmax*3);
ploty(1:3:end)=zeros(1,Nmax);
ploty(2:3:end)=da(ind(1:Nmax));
ploty(3:3:end)=nan;
plotx(1:3:end)=time(ind(1:Nmax));
plotx(2:3:end)=time(ind(1:Nmax));
plotx(3:3:end)=nan;
plot(plotx,ploty,'r-','linewidth',5);
xxlim=get(gca,'xlim');
plot(xxlim,[0 0],'k--');
plot(xxlim,[1 1],'m-.');
plot(xxlim,[-1 -1],'m-.');
yylim=get(gca,'ylim');
set(text(xxlim(2)/2,8/9*yylim(2),['$N_{CE,res}= ',num2str(NCE),'$']),'Interpreter','latex','fontsize',fontsize,'color','r');

plot([time(NCE) time(NCE)],yylim,'r--');

subplot(3,1,2);
plot(time,cumsum(de),'m-');hold all;
plot(etime,e,'r-');
plot(time,cumsum(di),'c-');
plot(etime,inc,'b-');
yylim=get(gca,'ylim');
plot([time(NCE) time(NCE)],yylim,'r--');
xlim(xxlim);

subplot(3,1,3);
plot(etime,phi,'k-');hold all;
ylim([0 360]);
yylim=get(gca,'ylim');
plot([time(NCE) time(NCE)],yylim,'r--');
xlim(xxlim);
