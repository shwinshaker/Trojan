function [aP,aT,eP,eT,IP,IT,RatioList,TimesList,MinDisList,SumList,RList,MaxCumsumList,EjectTime]=Fun_readMassHeap(PluName,DDir)

Dir='ServerMount'; 
di_name='di_record_inout';
de_name='de_record_inout';

RatioList1=1./exp(log(1):0.1:log(1000))';
RatioList1=sort(RatioList1);
RatioList2=exp(log(1):0.1:log(1000))';
RatioList=[RatioList1;RatioList2];

% RatioList1=exp(log(1e-3):0.1:log(1))';
% %RatioList1=sort(RatioList1);
% RatioList2=exp(log(1):0.1:log(1e3))';
% RatioList=[RatioList1;RatioList2];

% RatioList1=exp(log(1e-3):0.05:log(1))';
% %RatioList1=sort(RatioList1);
% RatioList2=exp(log(1):0.05:log(1e3))';
% RatioList=[RatioList1;RatioList2];

fname=cell(length(RatioList),1);
for i=1:length(RatioList1)
    Ratio=RatioList1(i);
    fname{i}=[PluName,'_',num2str(sprintf('%.4f',Ratio)),'MP'];
end
for i=1:length(RatioList2)
    Ratio=RatioList2(i);
    fname{length(RatioList1)+i}=[PluName,'_',num2str(sprintf('%.2f',Ratio)),'MP'];
end

MinDisList=zeros(length(RatioList),1);

aPL=zeros(length(RatioList),1);
IPL=zeros(length(RatioList),1);
ePL=zeros(length(RatioList),1);
eTL=zeros(length(RatioList),1);
aTL=zeros(length(RatioList),1);
ITL=zeros(length(RatioList),1);

TimesList=zeros(length(RatioList),1);

MaxList=zeros(length(RatioList),3);
MaxTime=zeros(length(RatioList),3);
SumList=zeros(length(RatioList),3);
StdList=zeros(length(RatioList),3);
RList=zeros(length(RatioList),3);
EjectTime=zeros(length(RatioList),1);
absDelList=zeros(length(RatioList),3);
MaxCumsumList=zeros(length(RatioList),3);

for i=1:length(RatioList)
    disp(fname{i});
    
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/plel.txt']);
    
    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if isempty(temp)
        ejectNo=size(tpel,1);
    else
        ejectNo=temp;
    end
    clear temp;
    ejecttime=tpel(ejectNo,1);
    
    tpel=tpel(1:ejectNo-1,:);
    plel=plel(1:ejectNo-1,:);
    
    aPL(i)=mean(plel(:,2));
    IPL(i)=mean(plel(:,4))/180*pi;
    ePL(i)=mean(plel(:,3));

    aTL(i)=mean(tpel(:,2));
    eTL(i)=mean(tpel(:,3));
    ITL(i)=mean(tpel(:,4))/180*pi;
    
    %absDelList(i,1)=abs(tpel(end,4)-tpel(1,4));
    absDelList(i,1)=max(abs(tpel(:,4)-mean(tpel(:,4))))/180*pi;
    %absDelList(i,2)=abs(tpel(end,3)-tpel(1,3));
    absDelList(i,2)=max(abs(tpel(:,3)-mean(tpel(:,3))));
    absDelList(i,3)=max(abs(tpel(:,2)-mean(tpel(:,2))));
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    
    di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',di_name,'.txt']);
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);
%     da_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',da_name,'.txt']);

    CEtime=CE_record(1:ejectCENo,1)/365.25;
    di_record_inout=di_record_inout(1:ejectCENo)/180*pi;
    de_record_inout=de_record_inout(1:ejectCENo);
%     da_record_inout=da_record_inout(1:ejectCENo);

    di_sum=cumsum(di_record_inout);
    de_sum=cumsum(de_record_inout);
%     da_sum=cumsum(da_record_inout);
    
    EjectTime(i)=ejecttime/365.25;
    TimesList(i)=length(di_record_inout);
    
    MinDisList(i)=min(CE_record(:,2));
    
    [MaxList(i,1),ind]=max(abs(di_record_inout));
    %MaxTime(i,1)=CEtime(abs(di_record_inout)==max(abs(di_record_inout)));
    MaxTime(i,1)=CEtime(ind);

    SumList(i,1)=abs(sum(di_record_inout));
    RList(i,1)=sum((di_record_inout-mean(di_record_inout)).^2);
    StdList(i,1)=std(di_record_inout);
    MaxCumsumList(i,1)=max(di_sum)-min(di_sum);
    
    [MaxList(i,2),ind]=max(abs(de_record_inout));
    MaxTime(i,2)=CEtime(ind);
    
    SumList(i,2)=abs(sum(de_record_inout));
    RList(i,2)=sum((de_record_inout-mean(de_record_inout)).^2);
    StdList(i,2)=std(de_record_inout);
    MaxCumsumList(i,2)=max(de_sum)-min(de_sum);
    
%     [MaxList(i,3),ind]=max(abs(da_record_inout));
%     MaxTime(i,3)=CEtime(ind);
%     
%     SumList(i,3)=abs(sum(da_record_inout));
%     RList(i,3)=sum((da_record_inout-mean(da_record_inout)).^2);
%     StdList(i,3)=std(da_record_inout);
%     MaxCumsumList(i,3)=max(da_sum)-min(da_sum);
end

aP=mean(aPL);
aT=mean(aTL);
eP=mean(ePL);
eT=mean(eTL);
IP=mean(IPL);
IT=mean(ITL);
 
