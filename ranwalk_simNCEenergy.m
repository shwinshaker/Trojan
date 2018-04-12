%% random simulation analysis
%% NCE and perturbating energy

if ~exist('TimesList','var') && ~exist('MinDisList','var')

clear;

Dir='swiftdata/Trojan/LAB/CE_realp/';
DDir='RanSimRanWalk2';
Name='1999CE119_2004UP10';

Nd=50;

MAXMIN=10;
TimesList=zeros(Nd,1);
SumDaList=zeros(Nd,1);

MinDisList=ones(Nd,MAXMIN);
MinTimesList=zeros(Nd,MAXMIN);

for i=1:Nd
    fname=[Name,'_da_',num2str(i)];
    disp(fname);
    tpel=load(['~/Documents/',Dir,DDir,'/',fname,'/tpel.txt']);
    
    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if isempty(temp)
        ejectNo=size(tpel,1);
    else
        ejectNo=temp;
    end
    clear temp;
    ejecttime=tpel(ejectNo,1);
    
    CE_record=load(['~/Documents/',Dir,DDir,'/',fname,'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    CE_record=CE_record(1:ejectCENo,:);
    
    di_record_inout=load(['~/Documents/',Dir,DDir,'/',fname,'/di_record_inout.txt']);
    da_record_inout=di_record_inout(1:ejectCENo,3);
    
    SumDaList(i)=abs(sum(da_record_inout));
    TimesList(i)=ejectCENo;
%     [MinDisList(i),MinTimesList(i)]=min(CE_record(:,2));
    
    Ind=find(CE_record(:,2)<0.0005);
    if length(Ind)>size(MinDisList,2)
        disp('Append Min List!');
        NcolNew=length(Ind)-size(MinDisList,2);
        Nrow=size(MinDisList,1);
        MinDisList=[MinDisList zeros(Nrow,NcolNew)];
        MinTimesList=[MinTimesList zeros(Nrow,NcolNew)];
    elseif isempty(Ind)
        disp('find No CE dis < 0.0005AU, pick the min dis!');
        [val,Ind]=min(CE_record(:,2));
    end
    MinDisList(i,1:length(Ind))=CE_record(Ind,2);
    MinTimesList(i,1:length(Ind))=Ind;
    
end

end

% shape=size(MinDisList);
% MinDisList(MinDisList==0)=1;
% MinDisList=reshape(MinDisList,shape);

figure;
set(gcf,'Position',[400,100,500,500],'color','w');
fontsize=15;
plot(0,0,'w+');hold all;
% plot(TimesList,min(MinDisList,[],2),'k.');
% plot(TimesList,SumDaList,'k.');
% plot(TimesList,MinTimesList,'k.');

[val,ind]=min(MinDisList,[],2);
MinMinTimesList=zeros(Nd,1);
for i=1:Nd
    MinMinTimesList(i)=MinTimesList(i,ind(i));
end
plot(TimesList,MinMinTimesList,'k.');

% plot(TimesList,max(MinTimesList,[],2),'r.','markersize',10);

% for i=1:Nd
%     plot(TimesList,MinTimesList(i,MinTimesList(i,:)>0 & MinTimesList(i,:)<max(MinTimesList(i,:),[],2)),'k.');
% end
xxlim=get(gca,'xlim');
yylim=get(gca,'ylim');
ulim=max(xxlim(2),yylim(2));
llim=min(xxlim(1),yylim(1));
plot([llim ulim],[llim ulim],'k-');
xlabel('$N_{CE,\,res}$','fontsize',fontsize,'interpreter','latex');
% ylabel('$A_{a,\,res}$','fontsize',fontsize,'interpreter','latex');
% ylabel('$D_{min}$','fontsize',fontsize,'interpreter','latex');
ylabel('$N_{min}$','fontsize',fontsize,'interpreter','latex');