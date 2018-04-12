%% total i change
clear;

ffname='RealPlutinos';
% fname1='1999CE119_1Gyr';
fname1='1999CE119_1Gyr';
fname2='2001KN77_1Gyr';
%fname1='1999CE119_1Gyr_10plumass';
%fname2='1998WV31_1Gyr';
%fname2='2001UO18_1Gyr';
%fname2='2001FU172_1Gyr';
%fname2='2001KN77_1Gyr_10plumass';
tag='di';

CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname1,'/CE_record.txt'));
di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname1,'/',tag,'_record_inout.txt'));

el=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname1,'/','tpel.txt'));

CE_record2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/CE_record.txt'));
di_record_inout2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/',tag,'_record_inout.txt'));
el2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname2,'/','tpel.txt'));

% CE_record2=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname2,'/CE_record.txt'));
% di_record_inout2=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname2,'/','di_record_inout.txt'));
% 
% lim2=0.2;
% lim1=0.4;

if strcmp(tag,'di')
    limdi=0.03;
    limsum=0.2;
    limDi=0.06;
    label='inc';
    unit=' /DEG';
    eldata=el(:,4);
    eldata2=el2(:,4);
else
    limdi=0.0005;
    limDi=0.002;
    limsum=0.01;
    label='e';
    unit='';
    eldata=el(:,3);
    eldata2=el2(:,3);
end

fontsize=15;

maxtime=max(max(el(:,1)),max(el2(:,1)));
%maxtime=max(el2(:,1));

di_sum=zeros(1,length(CE_record));
di_sum2=zeros(1,length(CE_record2));
di_sum_el=zeros(1,length(el)-1);
di_sum_el2=zeros(1,length(el2)-1);

for i=1:length(CE_record)
    di_sum(i)=sum(di_record_inout(1:i));
end

for i=1:length(el)-1
    di_sum_el(i)=el(i+1,4)-el(1,4);
end
for i=1:length(CE_record2)
    di_sum2(i)=sum(di_record_inout2(1:i));
end

for i=1:length(el2)-1
    di_sum_el2(i)=el2(i+1,4)-el2(1,4);
end
%% inc total change in each time interval
N=200;
Di=zeros(1,N);
Di2=zeros(1,N);
for i=1:N
    Di(i)=sum(di_record_inout((i-1)/N*maxtime<CE_record(:,1) & CE_record(:,1)<=i/N*maxtime));
end
for i=1:N
    Di2(i)=sum(di_record_inout2((i-1)/N*maxtime<CE_record2(:,1) & CE_record2(:,1)<=i/N*maxtime));
end


figure(1);
set(gcf,'Position',[400,100,900,600],'color','w');
subplot(4,1,1);
plot(el(:,1)/maxtime,eldata,'b-');hold on;
plot(el2(:,1)/maxtime,eldata2,'r-');
xlabel('Time','fontsize',fontsize);
ylabel([label,unit],'fontsize',fontsize);

subplot(4,1,2);
plot(CE_record(:,1)/maxtime,di_sum,'b-');hold on;
plot(CE_record2(:,1)/maxtime,di_sum2,'r-');hold on;
plot([0 1],[0 0],'k-');
ylim([-limsum limsum]);
xlabel('Time','fontsize',fontsize);
ylabel(['\delta',label,' sum',unit],'fontsize',fontsize);

subplot(4,1,3);
Ni=1:N;
Ni=Ni-1/2;
plot(Ni/N,Di,'b-');hold on;
plot(Ni/N,Di2,'r-');hold on;
plot([0 1],[0 0],'k-');

ylim([-limDi limDi]);
xlabel('Time','fontsize',fontsize);
ylabel(['\Delta',label,unit],'fontsize',fontsize);

subplot(4,1,4);
plot(CE_record(:,1)/maxtime,di_record_inout,'b.');hold on;
plot(CE_record2(:,1)/maxtime,di_record_inout2,'r.');hold on;
plot([0 1],[0 0],'k-');
ylim([-limdi limdi]);
xlabel('Time','fontsize',fontsize);
ylabel(['\delta',label,unit],'fontsize',fontsize);
%plot(el1(:,1)/maxtime,(el1(:,4)-el1(1,4))/max((el1(:,4)-el1(1,4))),'k-');hold on;
%plot(CE_record2(:,1)/maxtime,di_sum2,'r-');
