%% compare the tp inc in small step between fit&record
clear;
ffname='fit_record_contrast_201606';
% Dir_fit='1999CE119_1Gyr';
Dir_fit='1999CE119_tpAE_smlstp/fit';
%Dir_fit='1999CE119_forgedMass';
Dir_record='1999CE119_tpAE_smlstp';

datalines=254;

file_fit=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/AE_tp_smlstp.txt'];
fid_fit=fopen(file_fit,'r');
AE_fit=textscan(fid_fit,'%f %f %f %f %f %f','delimiter','\n');
fclose(fid_fit);
AE_fit=cell2mat(AE_fit);

dis_fit=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/dis_smlstp.txt'));

file_record=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/AE_tp_smlstp.txt'];
fid_record=fopen(file_record,'r');
textscan(fid_record,'%s',1,'delimiter','\n');
AE_record=textscan(fid_record,'%f %f %f %f %f %f',datalines,'delimiter','\n');
fclose(fid_record);
AE_record=cell2mat(AE_record);

file_record_dis=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/dis_smlstp.txt'];
fid_record_dis=fopen(file_record_dis,'r');
textscan(fid_record_dis,'%s',1,'delimiter','\n');
dis_record=textscan(fid_record_dis,'%f',datalines,'delimiter','\n');
fclose(fid_record_dis);
dis_record=cell2mat(dis_record);

%%pl el
file_fit_pl=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/AE_pl_smlstp.txt'];
fid_fit_pl=fopen(file_fit_pl,'r');
AE_fit_pl=textscan(fid_fit_pl,'%f %f %f %f %f %f','delimiter','\n');
fclose(fid_fit_pl);
AE_fit_pl=cell2mat(AE_fit_pl);

file_record_pl=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/AE_pl_smlstp.txt'];
fid_record_pl=fopen(file_record_pl,'r');
textscan(fid_record_pl,'%s',1,'delimiter','\n');
AE_record_pl=textscan(fid_record_pl,'%f %f %f %f %f %f',datalines,'delimiter','\n');
fclose(fid_record_pl);
AE_record_pl=cell2mat(AE_record_pl);

%% relative el
file_fit_rlt=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/AE_rlt_smlstp.txt'];
fid_fit_rlt=fopen(file_fit_rlt,'r');
AE_fit_rlt=textscan(fid_fit_rlt,'%f %f %f %f %f %f','delimiter','\n');
fclose(fid_fit_rlt);
AE_fit_rlt=cell2mat(AE_fit_rlt);

file_record_rlt=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/AE_rlt_smlstp.txt'];
fid_record_rlt=fopen(file_record_rlt,'r');
textscan(fid_record_rlt,'%s',1,'delimiter','\n');
AE_record_rlt=textscan(fid_record_rlt,'%f %f %f %f %f %f',datalines,'delimiter','\n');
fclose(fid_record_rlt);
AE_record_rlt=cell2mat(AE_record_rlt);


%%%%%%%%%%%%%%
r2hill=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/r2hill_record.txt'));
hill=sqrt(r2hill);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelFontSize=15;
figure(1);
set(gcf,'Position',[400,100,900,1000],'color','w');

subplot(3,3,1);
plot(AE_fit(:,2),AE_fit(:,3),'b.');hold all;axis square;grid on;
plot(AE_record(:,2),AE_record(:,3),'r.');
plot(AE_fit(ceil(length(AE_fit)/2),2),AE_fit(ceil(length(AE_fit)/2),3),'bo');
plot(AE_fit(ceil(length(AE_fit)/2)+1,2),AE_fit(ceil(length(AE_fit)/2)+1,3),'bo');
plot(AE_fit(1,2),AE_fit(1,3),'b*');

plot(AE_record(ceil(length(AE_record)/2),2),AE_record(ceil(length(AE_record)/2),3),'ro');
plot(AE_record(1,2),AE_record(1,3),'r*');
hold off;
xlabel('e','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);
legend({'Fit','Record'},'fontsize',15,'box','off');

subplot(3,3,2);
plot(AE_fit(:,1),AE_fit(:,3),'b.');hold all;axis square;grid on;
plot(AE_fit(ceil(length(AE_fit)/2),1),AE_fit(ceil(length(AE_fit)/2),3),'bo');
plot(AE_fit(ceil(length(AE_fit)/2)+1,1),AE_fit(ceil(length(AE_fit)/2)+1,3),'bo');
plot(AE_fit(1,1),AE_fit(1,3),'b*');

plot(AE_record(:,1),AE_record(:,3),'r.');
plot(AE_record(ceil(length(AE_record)/2),1),AE_record(ceil(length(AE_record)/2),3),'ro');
plot(AE_record(1,1),AE_record(1,3),'r*');
hold off;
xlabel('a /AU','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);

subplot(3,3,3);
plot(dis_fit/hill,AE_fit(:,3),'b.');hold all;axis square;grid on;
plot(dis_fit(1)/hill,AE_fit(1,3),'b*');

plot(dis_record/hill,AE_record(:,3),'r.');
plot(dis_record(1)/hill,AE_record(1,3),'r*');
hold off;
xlabel('CE distance /R_{hill}','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);

subplot(3,3,4);
plot(AE_fit_pl(:,2),AE_fit_pl(:,3),'b.');hold all;axis square;grid on;
plot(AE_record_pl(:,2),AE_record_pl(:,3),'r.');
plot(AE_fit_pl(ceil(length(AE_fit_pl)/2),2),AE_fit_pl(ceil(length(AE_fit_pl)/2),3),'bo');
plot(AE_fit_pl(ceil(length(AE_fit_pl)/2)+1,2),AE_fit_pl(ceil(length(AE_fit_pl)/2)+1,3),'bo');
plot(AE_fit_pl(1,2),AE_fit_pl(1,3),'b*');

plot(AE_record_pl(ceil(length(AE_record_pl)/2),2),AE_record_pl(ceil(length(AE_record_pl)/2),3),'ro');
plot(AE_record_pl(1,2),AE_record_pl(1,3),'r*');
hold off;
axis([0.2740 0.2742 1.7779 1.7781]);
xlabel('e','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);
legend({'Fit','Record'},'fontsize',15,'box','off');

subplot(3,3,5);
plot(AE_fit_pl(:,1),AE_fit_pl(:,3),'b.');hold all;axis square;grid on;
plot(AE_fit_pl(ceil(length(AE_fit_pl)/2),1),AE_fit_pl(ceil(length(AE_fit_pl)/2),3),'bo');
plot(AE_fit_pl(ceil(length(AE_fit_pl)/2)+1,1),AE_fit_pl(ceil(length(AE_fit_pl)/2)+1,3),'bo');
plot(AE_fit_pl(1,1),AE_fit_pl(1,3),'b*');

plot(AE_record_pl(:,1),AE_record_pl(:,3),'r.');
plot(AE_record_pl(ceil(length(AE_record_pl)/2),1),AE_record_pl(ceil(length(AE_record_pl)/2),3),'ro');
plot(AE_record_pl(1,1),AE_record_pl(1,3),'r*');
hold off;
axis([39.4162 39.4165 1.7779 1.7781]);
xlabel('a /AU','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);

subplot(3,3,7);
plot(AE_fit_rlt(:,2),AE_fit_rlt(:,3),'b.');hold all;axis square;grid on;
plot(AE_record_rlt(:,2),AE_record_rlt(:,3),'r.');
plot(AE_fit_rlt(ceil(length(AE_fit_rlt)/2),2),AE_fit_rlt(ceil(length(AE_fit_rlt)/2),3),'bo');
plot(AE_fit_rlt(ceil(length(AE_fit_rlt)/2)+1,2),AE_fit_rlt(ceil(length(AE_fit_rlt)/2)+1,3),'bo');
plot(AE_fit_rlt(1,2),AE_fit_rlt(1,3),'b*');

plot(AE_record_rlt(ceil(length(AE_record_rlt)/2),2),AE_record_rlt(ceil(length(AE_record_rlt)/2),3),'ro');
plot(AE_record_rlt(1,2),AE_record_rlt(1,3),'r*');
hold off;
xlabel('e','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);
legend({'Fit','Record'},'fontsize',15,'box','off');

subplot(3,3,8);
plot(AE_fit_rlt(:,1),AE_fit_rlt(:,3),'b.');hold all;axis square;grid on;
plot(AE_fit_rlt(ceil(length(AE_fit_rlt)/2),1),AE_fit_rlt(ceil(length(AE_fit_rlt)/2),3),'bo');
plot(AE_fit_rlt(ceil(length(AE_fit_rlt)/2)+1,1),AE_fit_rlt(ceil(length(AE_fit_rlt)/2)+1,3),'bo');
plot(AE_fit_rlt(1,1),AE_fit_rlt(1,3),'b*');

plot(AE_record_rlt(:,1),AE_record_rlt(:,3),'r.');
plot(AE_record_rlt(ceil(length(AE_record_rlt)/2),1),AE_record_rlt(ceil(length(AE_record_rlt)/2),3),'ro');
plot(AE_record_rlt(1,1),AE_record_rlt(1,3),'r*');
hold off;
xlabel('a /AU','fontsize',labelFontSize);
ylabel('inc /DEG','fontsize',labelFontSize);
