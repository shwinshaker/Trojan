%% random simulation analysis
%% NCE and perturbating energy

%% Input: CEOutput.txt

if ~exist('TimesList','var') && ~exist('MinDisList','var')

clear;

Dir='swiftdata/Trojan/LAB/CE_realp/';
DDir='RanSimRanWalk4';
Name='1999CE119_2004UP10';

Nd=1000;

TimesList=zeros(Nd,1);

SumList=zeros(Nd,6);
MaxList=zeros(Nd,6);

MinDisList=ones(Nd,1);
MinTimesList=zeros(Nd,1);

for i=1:Nd
    if i==9
        continue;
    end
    fname=[Name,'_da_',num2str(i)];
    disp(fname);
    fid=fopen(['~/Documents/',Dir,DDir,'/',fname,'/CEOutput.txt']);
    
    data=textscan(fid,'%d %f',1);
    TimesList(i)=data{1};
    
    data=textscan(fid,'%f %d %f',1);
    MinDisList(i)=data{1};
    MinTimesList(i)=data{2};
    
    data=textscan(fid,'%f',6);
    MaxList(i,:)=cell2mat(data)';
    data=textscan(fid,'%f',6);
    SumList(i,:)=cell2mat(data)';
    
    data=textscan(fid,'%f',1);
    
    fclose(fid); 
end

end

% shape=size(MinDisList);
% MinDisList(MinDisList==0)=1;
% MinDisList=reshape(MinDisList,shape);

%%% numerical example real particles

ffname='RealPlutinos';
fname={'1999CE119_1Gyr'};%'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};
MinDis1=zeros(length(fname),1);
Times1=zeros(length(fname),1);
for isub=1:length(fname)
    CE_record=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/CE_record.txt']);
    MinDis1(isub)=min(CE_record(:,2));
    Times1(isub)=length(CE_record(:,2));
end

ffname='RealPlutinosNpl';
fname={'1999CE119_1Gyr_5pl';'1999CE119_1Gyr_10pl';'1999CE119_1Gyr_15pl';'1999CE119_1Gyr_20pl';...,
    '1999CE119_1Gyr_25pl';'1999CE119_1Gyr_30pl';'1999CE119_1Gyr_35pl';...,
    '1999CE119_1Gyr_36pl';'1999CE119_1Gyr_37pl';'1999CE119_1Gyr_38pl';'1999CE119_1Gyr_39pl';'1999CE119_1Gyr_40pl'};
    %'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
MinDis2=zeros(length(fname),1);
Times2=zeros(length(fname),1);
for isub=1:length(fname)
    CE_record=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/CE_record.txt']);
    MinDis2(isub)=min(CE_record(:,2));
    Times2(isub)=length(CE_record(:,2));
end
MinDis=[MinDis1;MinDis2];
Times=[Times1;Times2];

figure;
set(gcf,'Position',[400,100,500,500],'color','w');

fontsize=15;

BottomRetainWidth=0.1;
LeftRetainWidth=0.12;
Height=0.42;
Width=0.8;

annotation('textbox',[Width/2-0.06 0.08+2*Height 0.05 0.05],'edgecolor','none','string',...,
    '1999CE119&2004UP10','fontweight','bold','fontsize',fontsize,'color','k');

annotation('textbox',[BottomRetainWidth+0.05 2*Height+0.02 0.05 0.05],'edgecolor','none','string',...,
    '$m_P = m_{Pluto}$','fontsize',fontsize,'color','r','Interpreter','latex');


loglog(0,0,'w+');hold all;
% plot(0,0,'w+');hold all;
hh3=loglog(TimesList,MinDisList,'k.');
% loglog(TimesList,abs(SumList(:,1)),'k.');
% [P,H]=polyfit(log(TimesList),log(abs(SumList(:,6))),1);
% loglog(TimesList,MinTimesList,'k.');

hh1=plot(Times1,MinDis1,'k.','markersize',20);
hh2=plot(Times2,MinDis2,'kx','markersize',8);


xxlim=get(gca,'xlim');
yylim=get(gca,'ylim');
loglog([mean(TimesList) mean(TimesList)],yylim,'k--','linewidth',2)
% ulim=max(xxlim(2),yylim(2));
% llim=min(xxlim(1),yylim(1));
% plot([llim ulim],[llim ulim],'k-');
aT=30;aP=40;gamma=3.5;mN=5.17e-5;mP=6.56e-9;
%FR=@(mu)((-2*aT^2*mu.^2+...,
%    2*(3^(1/3)*aP^2*aT^2*gamma^2*mu.^(8/3)+3*aT^4*mu.^4).^(1/2)/3^(1/2)).^(1/2));
% FR=@(mu)16/15*mu/(((2/aT-1/aP)^(1/2)-(1/aT)^(1/2))*(mN/aT)^(1/2));
%FR=@(mu)sqrt(2)*mu/(((2/aT-1/aP)^(1/2)-(1/aT)^(1/2))*(mN/aT)^(1/2));
FR=@(mu)Fun_Rk(aP,aT,mu/mP);
Rth=@(mu)(gamma*aP*(mu/3)^(1/3));
h1=plot(xxlim,[FR(mP) FR(mP)],'r-','linewidth',2);
xxx=exp(log(xxlim(1)):0.1:log(xxlim(2)));
h2=plot(xxx,Rth(mP)*(1./xxx).^(1/2),'b-','linewidth',2);

legend([h1 h2],{'$R_k$','$R_\varepsilon$'},'fontsize',fontsize,'Interpreter','latex','location','southwest');

xlabel('$N_{CE,\,res}$','fontsize',fontsize,'interpreter','latex');
% ylabel('$A_{a,\,res}$','fontsize',fontsize,'interpreter','latex');
ylabel('$D_{min}\rm(AU)$','fontsize',fontsize,'interpreter','latex');
% ylabel('$N_{min}$','fontsize',fontsize,'interpreter','latex');

ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,[hh1 hh2 hh3],{'NS','ClonedNS','Mixed'},'fontsize',fontsize,...,
'location','northeast');
