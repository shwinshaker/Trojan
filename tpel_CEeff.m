%tpel_CEeff
clear;
ffname='RealTrojans';
tp_name={'2004UP10';'2005TN53';'2005TO74';'2006RJ103';'2007VL305';...,
    '2004KV18';'2008LC18';'2010TS191';'2010TT191';'2011HM102';'2011SO277';...,
    '2011WG157';'2012UV177';'2013KY18';'2014QO441'};
CE_times=cell(length(tp_name),1);
inc=cell(length(tp_name),1);
%amp=cell(length(tp_name),1);
Di=cell(length(tp_name),1);
Dir='ServerMount'; % 'swiftdata'
for i=1:length(tp_name)
    fname=[tp_name{i},'_1Gyr'];
    disp(tp_name{i});
    %CE_record=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
    tpElIn=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',tp_name{i},'.txt'));
    tp_el=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt'));
    di_record_inout=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/','di_record_inout.txt'));
    %eval(['CE_record',num2str(i),'=CE_record;']);
    %eval(['CE_time',num2str(i),'=length(CE_record);']);
    %eval(['pl_el',num2str(i),'=pl_el;']);
    CE_times{i}=length(di_record_inout);
    inc{i}=mean(tp_el(:,4));
    %amp{i}=tpElIn(9);
    Di{i}=sum(di_record_inout);
end
datacell=[tp_name CE_times inc Di];%amp
datacell=sortrows(datacell,-2);
maxDi=max(abs(cell2mat(datacell(:,4))));
fontsize=15;
data=cell2mat(datacell(:,2:end));

%%% Trojan info

filename='~/Desktop/NTData1.txt';
disp(['Trojan info file: ',filename]);
fid=fopen(filename,'r');
headerlines=textscan(fid,'%s',1,'delimiter','\n');
info=textscan(fid,'%s %f %f %f %f %f %f %s %*[^\n]');
fclose(fid);

figure(1);
set(gcf,'Position',[400,100,800,500],'color','w');

subplot(1,2,1);
plot(0,0,'w');hold all;
for i=1:length(data)
    infoindex=find(strcmp(info{1},datacell{i,1}));
    if strcmp(info{8}{infoindex},'L4')
        rp=plot(data(i,2),data(i,1),'r+');
    elseif strcmp(info{8}{infoindex},'L5')
        bp=plot(data(i,2),data(i,1),'b+');
    else
        disp('error!');
    end
end
hold off;
xlabel('$$\rm{\overline{Inc.}}$$','Interpreter','latex','fontsize',fontsize);
ylabel('N CE','fontsize',fontsize);
legend([rp(1),bp(1)],{'L4','L5'},'fontsize',fontsize);

subplot(1,2,2);
plot(0,0,'w');hold all;
for i=1:length(data)
    infoindex=find(strcmp(info{1},datacell{i,1}));
    if strcmp(info{8}{infoindex},'L4')
        rp=plot(data(i,2),abs(data(i,3)),'r+');
    elseif strcmp(info{8}{infoindex},'L5')
        bp=plot(data(i,2),abs(data(i,3)),'b+');
    else
        disp('error!');
    end
end
xlabel('$$\rm{\overline{Inc.}}$$','Interpreter','latex','fontsize',fontsize);
ylabel('|\Delta Inc.|','fontsize',fontsize);
legend([rp(1),bp(1)],{'L4','L5'},'fontsize',fontsize);
%legend('boxoff');
% subplot(2,2,3);
% plot(data(:,3),data(:,1),'k+');
% subplot(2,2,4);
% plot(data(:,3),abs(data(:,4)),'k+');