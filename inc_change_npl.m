%% total i change
clear;

%ffname='RealPlutinosNpl';
ffname='MassTest2';
fname='1999CE119_12.18MP';
tag='di';

% CE_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
% di_record_inout=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','di_record_perturb.txt'));
% el=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','tpel.txt'));

CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
%CE_record(54,:)=[]; %% CO-CE
ddata=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt'));
el=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/','tpel.txt'));


%ddata(isnan(ddata))=0;
if strcmp(tag,'di')
%     limdi=0.5;%0.03;
%     limsum=1.0;%0.2;
%     limDi=0.5;
    limdi=0.1;
    limsum=0.2;
    limDi=0.1;
    label='inc';
    unit=' /DEG';
    eldata=el(:,4);
else
    limdi=0.001;
    limDi=0.001;
    limsum=0.02;
    label='e';
    unit='';
    eldata=el(:,3);
end

fontsize=15;
N=200;

eltime=el(:,1);
maxtime=max(el(:,1));
time=CE_record(:,1); %% time scale

di_sum=cumsum(ddata);

%% inc total change in each time interval

Di=zeros(1,N);
for i=1:N
    Di(i)=sum(ddata((i-1)/N*maxtime<time & time<=i/N*maxtime));
end

figure(2);
set(gcf,'Position',[400,100,900,600],'color','w');

row=5;

subplot(row,1,1);

% [elpeak,elpeaktime]=findpeaks(eldata,eltime);
% [eltrough,eltroughtime]=findpeaks(-eldata,eltime);
% eltrough=-eltrough;
% 
% %elsplinetime=[elpeaktime(2:end)-diff(elpeaktime)/4;elpeaktime(2:end)-diff(elpeaktime)*3/4];
% if length(elpeaktime)==length(eltroughtime)-1
%     elsplinetime=[(elpeaktime+eltroughtime(2:end))/2;(elpeaktime+eltroughtime(1:end-1))/2];
% elseif length(elpeaktime)==length(eltroughtime)+1
%     elsplinetime=[(elpeaktime(2:end)+eltroughtime)/2;(elpeaktime(1:end-1)+eltroughtime)/2];
% else
%     quit;
% end
% %elsplinetime=[(elpeaktime+eltroughtime(2:end))/2;(elpeaktime+eltroughtime(1:end-1))/2];
% elsplinetime=sort(elsplinetime);
% 
% rank=20;
% hpoly=polyfit(eltime,eldata,rank);
% elspline=polyval(hpoly,elsplinetime);
% % elspline=spline(eltime,eldata,elsplinetime);

plot(eltime/maxtime,el(:,2),'r-');
% plot(eltime/maxtime,eldata,'b-');hold all;
% plot(elpeaktime/maxtime,elpeak,'b--');
% plot(eltroughtime/maxtime,eltrough,'b--');
% plot(elsplinetime/maxtime,elspline,'r-');
hold off;
xlabel('Time','fontsize',fontsize);
ylabel([label,unit],'fontsize',fontsize);

subplot(row,1,2);

% hpoly=polyfit(time,di_sum,rank);
% displine=polyval(hpoly,elsplinetime);

% displine=spline(time,di_sum,elsplinetime);

% plot(time/maxtime,di_sum,'b-');hold on;
plot(eltime/maxtime,el(:,3),'r-');hold on;
% plot(elsplinetime/maxtime,displine,'r-');
plot([0 1],[0 0],'k-');

% ylim([-limsum limsum]);
xlabel('Time','fontsize',fontsize);
ylabel('e','fontsize',fontsize);


subplot(row,1,3);
plot(time/maxtime,di_sum,'b-');hold on;
plot([0 1],[0 0],'k-');

% ylim([-limsum limsum]);
xlabel('Time','fontsize',fontsize);
ylabel(['\delta',label,' sum',unit],'fontsize',fontsize);


subplot(row,1,4); %% bar is better?
Ni=1:N;
Ni=Ni-1/2;
% plot(Ni/N,Di,'b-');hold on;
bar(Ni/N,Di,'b');hold on;
plot([0 1],[0 0],'k-');

% ylim([-limDi limDi]);
xlabel('Time','fontsize',fontsize);
ylabel(['\Delta',label,unit],'fontsize',fontsize);

subplot(row,1,5);
plot(time/maxtime,ddata,'b.');hold on;
plot([0 1],[0 0],'k-');
% ylim([-limdi limdi]);
xlabel('Time','fontsize',fontsize);
ylabel(['\delta',label,unit],'fontsize',fontsize);
%plot(el1(:,1)/maxtime,(el1(:,4)-el1(1,4))/max((el1(:,4)-el1(1,4))),'k-');hold on;
%plot(CE_record2(:,1)/maxtime,di_sum2,'r-');
