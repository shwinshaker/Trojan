%% fit all

[PT1,HT1]=polyfit(log(RatioList1),log(TimesList(1:length(RatioList1),1)),1);  
RT1=corrcoef(log(RatioList1),log(TimesList(1:length(RatioList1),1)));
Timesfit1=exp(polyval(PT1,log(RatioList1)));

[PT2,HT2]=polyfit(log(RatioList2),log(TimesList(length(RatioList1)+1:end,1)),1);  
RT2=corrcoef(log(RatioList2),log(TimesList(length(RatioList1)+1:end,1)));
Timesfit2=exp(polyval(PT2,log(RatioList2)));

[PSum1,HSum1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),1)),1);  
RSum1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),1)));
Sumfit1=exp(polyval(PSum1,log(RatioList1)));

[PSum2,HSum2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,1)),1);  
RSum2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,1)));
Sumfit2=exp(polyval(PSum2,log(RatioList2)));

%%% Cum Sum Max
[PCumsum1,HCumsum1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)),1);  
RCumsum1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)));
Cumsumfit1=exp(polyval(PCumsum1,log(RatioList1)));

[PCumsum2,HCumsum2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)),1);  
RCumsum2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)));
Cumsumfit2=exp(polyval(PCumsum2,log(RatioList2)));

%%%
[PStd1,HStd1]=polyfit(log(RatioList1),log(StdList(1:length(RatioList1),1)),1);  
RStd1=corrcoef(log(RatioList1),log(StdList(1:length(RatioList1),1)));
Stdfit1=exp(polyval(PStd1,log(RatioList1)));

[PStd2,HStd2]=polyfit(log(RatioList2),log(StdList(length(RatioList1)+1:end,1)),1);  
RStd2=corrcoef(log(RatioList2),log(StdList(length(RatioList1)+1:end,1)));
Stdfit2=exp(polyval(PStd2,log(RatioList2)));

[PR1,HR1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),1)),1);  
RR1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),1)));
Rfit1=exp(polyval(PR1,log(RatioList1)));

[PR2,HR2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,1)),1);  
RR2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,1)));
Rfit2=exp(polyval(PR2,log(RatioList2)));

%%de
[PSumde1,HSumde1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),2)),1);  
RSumde1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),2)));
Sumfitde1=exp(polyval(PSumde1,log(RatioList1)));

[PSumde2,HSumde2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,2)),1);  
RSumde2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,2)));
Sumfitde2=exp(polyval(PSumde2,log(RatioList2)));

%% Cum Sum Max
[PCumsumde1,HCumsumde1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)),1);  
RCumsumde1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)));
Cumsumfitde1=exp(polyval(PCumsumde1,log(RatioList1)));

[PCumsumde2,HCumsumde2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)),1);  
RCumsumde2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)));
Cumsumfitde2=exp(polyval(PCumsumde2,log(RatioList2)));


[PRde1,HRde1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),2)),1);  
RRde1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),2)));
Rfitde1=exp(polyval(PRde1,log(RatioList1)));

[PRde2,HRde2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,2)),1);  
RRde2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,2)));
Rfitde2=exp(polyval(PRde2,log(RatioList2)));

%%da
[PSumda1,HSumda1]=polyfit(log(RatioList1),log(SumList(1:length(RatioList1),3)),1);  
RSumda1=corrcoef(log(RatioList1),log(SumList(1:length(RatioList1),3)));
Sumfitda1=exp(polyval(PSumda1,log(RatioList1)));

[PSumda2,HSumda2]=polyfit(log(RatioList2),log(SumList(length(RatioList1)+1:end,3)),1);  
RSumda2=corrcoef(log(RatioList2),log(SumList(length(RatioList1)+1:end,3)));
Sumfitda2=exp(polyval(PSumda2,log(RatioList2)));

%% Cum Sum Max
[PCumsumda1,HCumsumda1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),3)),1);  
RCumsumda1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),3)));
Cumsumfitda1=exp(polyval(PCumsumda1,log(RatioList1)));

[PCumsumda2,HCumsumda2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,3)),1);  
RCumsumda2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,3)));
Cumsumfitda2=exp(polyval(PCumsumda2,log(RatioList2)));


[PRda1,HRda1]=polyfit(log(RatioList1),log(RList(1:length(RatioList1),3)),1);  
RRda1=corrcoef(log(RatioList1),log(RList(1:length(RatioList1),3)));
Rfitda1=exp(polyval(PRda1,log(RatioList1)));

[PRda2,HRda2]=polyfit(log(RatioList2),log(RList(length(RatioList1)+1:end,3)),1);  
RRda2=corrcoef(log(RatioList2),log(RList(length(RatioList1)+1:end,3)));
Rfitda2=exp(polyval(PRda2,log(RatioList2)));
