%% total i change
clear;

ffname='fit_record_contrast_201606';
% fname1='1999CE119_1Gyr';
fname='1999CE119_1Gyr';
% fname2='2001KN77_1Gyr';


CE_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
di_record_inout=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','di_record_inout.txt'));

el=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','tpel.txt'));

ran_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/ran_record.txt'));
di_ran_fit_ptb=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','ran_di_fit_ptb.txt'));

% CE_record2=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname2,'/CE_record.txt'));
% di_record_inout2=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname2,'/','di_record_inout.txt'));
% 
lim2=0.01;
lim1=0.2;
fontsize=15;

maxtime=max(el(:,1));
ran_time=ran_record(:,1)/ran_record(end,1)*maxtime; %% time scale
% maxtime=max(max(CE_record1(:,1)),max(CE_record2(:,1)));
di_sum=zeros(1,length(CE_record));
di_sum_ranfit=zeros(1,length(ran_record));
% di_sum2=zeros(1,length(CE_record2));
di_sum_el=zeros(1,length(el)-1);

for i=1:length(CE_record)
    di_sum(i)=sum(di_record_inout(1:i));
end

for i=1:length(ran_record)
    di_sum_ranfit(i)=sum(di_ran_fit_ptb(1:i));
end

for i=1:length(el)-1
    di_sum_el(i)=el(i+1,4)-el(1,4);
end
% for i=1:length(CE_record2)
%     di_sum2(i)=sum(di_record_inout2(1:i));
% end
%% inc total change in each time interval
N=200;
Di=zeros(1,N);
Diran=zeros(1,N);
for i=1:N
    Di(i)=sum(di_record_inout((i-1)/N*maxtime<CE_record(:,1) & CE_record(:,1)<=i/N*maxtime));
end
for i=1:N
    Diran(i)=sum(di_ran_fit_ptb((i-1)/N*maxtime<ran_time & ran_time<=i/N*maxtime));
end

figure(1);
set(gcf,'Position',[400,100,900,600],'color','w');
subplot(4,1,1);
plot(el(2:end,1)/maxtime,di_sum_el,'k-');
xlabel('Time','fontsize',fontsize);
ylabel('\deltainc hitherto','fontsize',fontsize);

subplot(4,1,2);
plot(CE_record(:,1)/maxtime,di_sum,'b-');hold on;
plot(ran_time/maxtime,di_sum_ranfit,'r-');
ylim([-lim1 lim1]);
xlabel('Time','fontsize',fontsize);
ylabel('\deltainc sum','fontsize',fontsize);

subplot(4,1,3);
Ni=1:N;
Ni=Ni-1/2;
plot(Ni/N,Di,'b-');hold on;
plot(Ni/N,Diran,'r-');hold on;
xlabel('Time','fontsize',fontsize);
ylabel('\Deltainc','fontsize',fontsize);

subplot(4,1,4);
plot(CE_record(:,1)/maxtime,di_record_inout,'b.');hold on;
plot(ran_time/maxtime,di_ran_fit_ptb,'r.');
plot([0 1],[0 0],'k-');
ylim([-lim2 lim2]);
xlabel('Time','fontsize',fontsize);
ylabel('\deltainc','fontsize',fontsize);
%plot(el1(:,1)/maxtime,(el1(:,4)-el1(1,4))/max((el1(:,4)-el1(1,4))),'k-');hold on;
%plot(CE_record2(:,1)/maxtime,di_sum2,'r-');
