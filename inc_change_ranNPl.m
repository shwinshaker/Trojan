%% total i change
clear;

ffname='RanPlutinos';
% fname1='1999CE119_1Gyr';
fname='1999CE119_500k';
%fname='2001KN77_1Gyr';

ran_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ran_record.txt'));
di_ran_fit_ptb=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/','di_fit_perturb.txt'));
ierror=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/','ierror.txt'));
ran_record(ierror,:)=[];

lim2=0.1;
lim1=1.0;
fontsize=15;

maxtime=3.65e11;
ran_time=ran_record(:,1)/ran_record(end,1)*maxtime; %% time scale

% for i=1:length(ran_record)
%     di_sum_ranfit(i)=sum(di_ran_fit_ptb(1:i));
% end
di_sum_ranfit=cumsum(di_ran_fit_ptb);

%% inc total change in each time interval
N=50;

Diran=zeros(N,1);
for i=1:N
    Diran(i)=sum(di_ran_fit_ptb((i-1)/N*maxtime<ran_time & ran_time<=i/N*maxtime));
end

figure(1);
set(gcf,'Position',[400,100,900,600],'color','w');

subplot(3,1,1);
plot(ran_time/maxtime,di_sum_ranfit,'r-');
ylim([-lim1 lim1]);
xlabel('Time','fontsize',fontsize);
ylabel('\deltainc sum','fontsize',fontsize);

subplot(3,1,2);
Ni=1:N;
Ni=Ni-1/2;
plot(Ni/N,Diran,'r.-');hold on;
xlabel('Time','fontsize',fontsize);
ylabel('\Deltainc','fontsize',fontsize);

subplot(3,1,3);
plot(ran_time/maxtime,di_ran_fit_ptb,'r.');hold on;
plot([0 1],[0 0],'k-');
ylim([-lim2 lim2]);
xlabel('Time','fontsize',fontsize);
ylabel('\deltainc','fontsize',fontsize);
%plot(el1(:,1)/maxtime,(el1(:,4)-el1(1,4))/max((el1(:,4)-el1(1,4))),'k-');hold on;
%plot(CE_record2(:,1)/maxtime,di_sum2,'r-');
