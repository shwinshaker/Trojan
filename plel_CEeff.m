%% plutino initial el ---- CE effects (better for average el)
%% 2003WU172 run away
clear;
ffname='RealPlutinos';
pl_name={'2001KN77';'1999CE119';'1998WV31';'2001FU172';'1993SB'; ...,
    '2001UO18';'1994TB';'1999CM158';'2000YH2';'2002CE251';'2000FB8';'1996TP66';'1998HQ151';...,
    '2001RX143';'1995QY9';'1996SZ4';...,
    '2000GE147';'2001KD77';'2001FR185';'1998HK151';...,
    '2002CW224';'2001KY76';'2000CK105';'1995HM5'};
CE_times=cell(length(pl_name),1);
ecc=cell(length(pl_name),1);
inc=cell(length(pl_name),1);
amp=cell(length(pl_name),1);
Di=cell(length(pl_name),1);
Dir='ServerMount'; % 'swiftdata'
for i=1:length(pl_name)
    fname=[pl_name{i},'_1Gyr'];
    disp(pl_name{i});
    CE_record=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
    plElIn=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',pl_name{i},'.txt'));
    pl_el=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/plel.txt'));
    di_record_inout=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/','di_record_inout.txt'));
    %eval(['CE_record',num2str(i),'=CE_record;']);
    %eval(['CE_time',num2str(i),'=length(CE_record);']);
    %eval(['pl_el',num2str(i),'=pl_el;']);
    CE_times{i}=length(CE_record);
    ecc{i}=mean(pl_el(:,3));
    inc{i}=mean(pl_el(:,4));
    amp{i}=plElIn(9);
    Di{i}=sum(di_record_inout);
end
datacell=[pl_name CE_times ecc inc amp Di];
datacell=sortrows(datacell,-2);
maxDi=max(abs(cell2mat(datacell(:,6))));
fontsize=15;

basicSize=150;
maxtimes=3000;
markersize=zeros(length(pl_name),1);
for i=1:length(pl_name)
    markersize(i)=basicSize/(1-log(abs(datacell{i,6})/maxDi));
end

figure(1);
set(gcf,'Position',[400,100,1000,500],'color','w');

subplot(1,2,1);
plot(0,0);hold all;
for i=1:length(pl_name)
%     eval(['x=pl_el',num2str(i),'(2);']);
%     eval(['y=pl_el',num2str(i),'(3);']);
%     eval(['CE_time=CE_time',num2str(i),';']);
%     plot(x,y,'k.','markersize',CE_time/maxtimes*basicSize);
    if markersize(i)~=0
        plot(datacell{i,3},datacell{i,4},'.','markersize',markersize(i),'color',[datacell{i,2}/maxtimes 0 1-datacell{i,2}/maxtimes]);hold on;
    else
        plot(datacell{i,3},datacell{i,4},'.','color',[datacell{i,2}/maxtimes 0 1-datacell{i,2}/maxtimes]);hold on;
    end

end
hold off;
xlabel('$$\rm{\overline{e}}$$','Interpreter','latex','fontsize',fontsize);
ylabel('$$\rm{\overline{i}\ /DEG}$$','Interpreter','latex','fontsize',fontsize);

subplot(1,2,2);
plot(0,0);hold all;
for i=1:length(pl_name)
%     eval(['x=pl_el',num2str(i),'(2);']);
%     eval(['y=pl_el',num2str(i),'(3);']);
%     eval(['CE_time=CE_time',num2str(i),';']);
%     plot(x,y,'k.','markersize',CE_time/maxtimes*basicSize);
    if markersize(i)~=0
        plot(datacell{i,5},datacell{i,4},'.','markersize',markersize(i),'color',[datacell{i,2}/maxtimes 0 1-datacell{i,2}/maxtimes]);hold on;
    else
        plot(datacell{i,5},datacell{i,4},'.','color',[datacell{i,2}/maxtimes 0 1-datacell{i,2}/maxtimes]);hold on;
    end
end
hold off;
xlim([0 180]);
xlabel('$$\rm{\overline{\Delta \phi}\ /DEG}$$','Interpreter','latex','fontsize',fontsize);
ylabel('$$\rm{\overline{i}\ /DEG}$$','Interpreter','latex','fontsize',fontsize);