%% small step xy position
clear;
ffname='fit_record_contrast_201606';
Dir_fit='1999CE119_tpAE_smlstp/fit';
Dir_record='1999CE119_tpAE_smlstp';

file_fit_tp=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/XV_tp_smlstp.txt'];
fid_fit_tp=fopen(file_fit_tp,'r');
XV_fit_tp=textscan(fid_fit_tp,'%f %f %f %f %f %f','delimiter','\n');
fclose(fid_fit_tp);
XV_fit_tp=cell2mat(XV_fit_tp);

file_fit_pl=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_fit,'/XV_pl_smlstp.txt'];
fid_fit_pl=fopen(file_fit_pl,'r');
XV_fit_pl=textscan(fid_fit_pl,'%f %f %f %f %f %f','delimiter','\n');
fclose(fid_fit_pl);
XV_fit_pl=cell2mat(XV_fit_pl);

datalines=254;
file_record_tp=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/XV_tp_smlstp.txt'];
fid_record_tp=fopen(file_record_tp,'r');
textscan(fid_record_tp,'%s',1,'delimiter','\n');
XV_record_tp=textscan(fid_record_tp,'%f %f %f %f %f %f',datalines,'delimiter','\n');
fclose(fid_record_tp);
XV_record_tp=cell2mat(XV_record_tp);

file_record_pl=['~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',Dir_record,'/XV_pl_smlstp.txt'];
fid_record_pl=fopen(file_record_pl,'r');
textscan(fid_record_pl,'%s',1,'delimiter','\n');
XV_record_pl=textscan(fid_record_pl,'%f %f %f %f %f %f',datalines,'delimiter','\n');
fclose(fid_record_pl);
XV_record_pl=cell2mat(XV_record_pl);

figure(1);
set(gcf,'Position',[400,100,1000,500],'color','w');

% ????tp pl?? ????
subplot(1,2,1);
%% tp
plot(XV_fit_tp(:,1),XV_fit_tp(:,2),'b.');hold all;axis square;grid on;
plot(XV_record_tp(:,1),XV_record_tp(:,2),'r.');
%% pl
plot(XV_fit_pl(:,1),XV_fit_pl(:,2),'b--');
plot(XV_record_pl(:,1),XV_record_pl(:,2),'r--');
legend({'Trojan-fit','Trojan-record','Plutino-fit','Plutino-record'},'fontsize',15,'box','off');

%%% tp
plot(XV_fit_tp(ceil(length(XV_fit_tp)/2),1),XV_fit_tp(ceil(length(XV_fit_tp)/2),2),'bo');
plot(XV_fit_tp(ceil(length(XV_fit_tp)/2)+1,1),XV_fit_tp(ceil(length(XV_fit_tp)/2)+1,2),'bo');
plot(XV_fit_tp(1,1),XV_fit_tp(1,2),'b*');

plot(XV_record_tp(ceil(length(XV_record_tp)/2),1),XV_record_tp(ceil(length(XV_record_tp)/2),2),'ro');
plot(XV_record_tp(1,1),XV_record_tp(1,2),'r*');

%%% pl 

plot(XV_fit_pl(ceil(length(XV_fit_pl)/2),1),XV_fit_pl(ceil(length(XV_fit_pl)/2),2),'bo');
plot(XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,1),XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,2),'bo');
plot(XV_fit_pl(1,1),XV_fit_pl(1,2),'b*');

plot(XV_record_pl(ceil(length(XV_record_pl)/2),1),XV_record_pl(ceil(length(XV_record_pl)/2),2),'ro');
plot(XV_record_pl(1,1),XV_record_pl(1,2),'r*');
hold off;
xlabel('x /AU','fontsize',15);
ylabel('y /AU','fontsize',15);


subplot(1,2,2);
%% relative
plot(0,0,'b+');hold all;axis square;grid on;
plot(XV_fit_tp(:,1)-XV_fit_pl(:,1),XV_fit_tp(:,2)-XV_fit_pl(:,2),'b-');
plot(XV_fit_tp(ceil(length(XV_fit_pl)/2),1)-XV_fit_pl(ceil(length(XV_fit_pl)/2),1),XV_fit_tp(ceil(length(XV_fit_pl)/2),2)-XV_fit_pl(ceil(length(XV_fit_pl)/2),2),'bo');
plot(XV_fit_tp(ceil(length(XV_fit_pl)/2)+1,1)-XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,1),XV_fit_tp(ceil(length(XV_fit_pl)/2)+1,2)-XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,2),'bo');
plot(XV_fit_tp(1,1)-XV_fit_pl(1,1),XV_fit_tp(1,2)-XV_fit_pl(1,2),'b*');

plot(XV_record_tp(:,1)-XV_record_pl(:,1),XV_record_tp(:,2)-XV_record_pl(:,2),'r-');
plot(XV_record_tp(ceil(length(XV_record_pl)/2),1)-XV_record_pl(ceil(length(XV_record_pl)/2),1),XV_record_tp(ceil(length(XV_record_pl)/2),2)-XV_record_pl(ceil(length(XV_record_pl)/2),2),'ro');
plot(XV_record_tp(ceil(length(XV_record_pl)/2)+1,1)-XV_record_pl(ceil(length(XV_record_pl)/2)+1,1),XV_record_tp(ceil(length(XV_record_pl)/2)+1,2)-XV_record_pl(ceil(length(XV_record_pl)/2)+1,2),'ro');
plot(XV_record_tp(1,1)-XV_record_pl(1,1),XV_record_tp(1,2)-XV_record_pl(1,2),'r*');hold off;
xlabel('x /AU','fontsize',15);
ylabel('y /AU','fontsize',15);

% subplot(2,2,2);
% plot(0,0,'b+');hold all;axis square;grid on;
% plot(XV_fit_tp(:,1)-XV_fit_pl(:,1),XV_fit_tp(:,2)-XV_fit_pl(:,2),'b-');
% plot(XV_fit_tp(ceil(length(XV_fit_pl)/2),1)-XV_fit_pl(ceil(length(XV_fit_pl)/2),1),XV_fit_tp(ceil(length(XV_fit_pl)/2),2)-XV_fit_pl(ceil(length(XV_fit_pl)/2),2),'bo');
% plot(XV_fit_tp(ceil(length(XV_fit_pl)/2)+1,1)-XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,1),XV_fit_tp(ceil(length(XV_fit_pl)/2)+1,2)-XV_fit_pl(ceil(length(XV_fit_pl)/2)+1,2),'bo');
% plot(XV_fit_tp(1,1)-XV_fit_pl(1,1),XV_fit_tp(1,2)-XV_fit_pl(1,2),'b*');
% 
% plot(XV_record_tp(:,1)-XV_record_pl(:,1),XV_record_tp(:,2)-XV_record_pl(:,2),'r-');
% plot(XV_record_tp(ceil(length(XV_record_pl)/2),1)-XV_record_pl(ceil(length(XV_record_pl)/2),1),XV_record_tp(ceil(length(XV_record_pl)/2),2)-XV_record_pl(ceil(length(XV_record_pl)/2),2),'ro');
% plot(XV_record_tp(ceil(length(XV_record_pl)/2)+1,1)-XV_record_pl(ceil(length(XV_record_pl)/2)+1,1),XV_record_tp(ceil(length(XV_record_pl)/2)+1,2)-XV_record_pl(ceil(length(XV_record_pl)/2)+1,2),'ro');
% plot(XV_record_tp(1,1)-XV_record_pl(1,1),XV_record_tp(1,2)-XV_record_pl(1,2),'r*');hold off;
