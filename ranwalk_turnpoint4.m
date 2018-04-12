%%MassTest
if exist('TimesList','var') 
    tag='repeat';
else
    tag='data';
end

%if ~strcmp(tag,'repeat')
clear;
Dir='ServerMount'; 

di_name='di_record_inout';
de_name='de_record_inout';

figure;
BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.4;

set(gcf,'Position',[400,100,700,500],'color','w');
row=2;
col=2;

linewidth=2;
fontsize=15;

xxlim=[1e-3 1e3];
xxtick=power(10,-3:1:3);
xp=1.5e-3;

ReviseLineStyle='-';
ReviseColor='c';


for ipl=1:4
    switch ipl
        case 1
            PluName='1999CE119';
            DDir='MassTest2';
            titlename='1999CE119&2004UP10';
        case 2
            PluName='2001FU172';
            DDir='MassTest3';
            titlename='2001FU172&2004UP10';
        case 3
            PluName='1999CE119_2006RJ103';
            DDir='MassTest4';
            titlename='1999CE119&2006RJ103';
        case 4
            PluName='2001FU172_2006RJ103';
            DDir='MassTest5';
            titlename='2001FU172&2006RJ103';
    end

    %%%%%%%%%%%%% Revise 1&2 %%%%%%%%%%%%%
    if strcmp(PluName,'1999CE119') 
        barrier=5.0;
    elseif strcmp(PluName,'1999CE119_2006RJ103') 
        barrier=10.0;
    else
        barrier=30.0;
    end
    
    RatioList1=1./exp(log(1):0.1:log(1000))';
    RatioList1=sort(RatioList1);
    RatioList2=exp(log(1):0.1:log(1000))';
    RatioList=[RatioList1;RatioList2];

    fname=cell(length(RatioList),1);
    for i=1:length(RatioList1)
        Ratio=RatioList1(i);
        fname{i}=[PluName,'_',num2str(sprintf('%.4f',Ratio)),'MP'];
    end
    for i=1:length(RatioList2)
        Ratio=RatioList2(i);
        fname{length(RatioList1)+i}=[PluName,'_',num2str(sprintf('%.2f',Ratio)),'MP'];
    end
   
    iF=find(RatioList2<=barrier,1,'last');
    RatioList1=[RatioList1;RatioList2(1:iF)];
    RatioList2(1:iF)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %EjectTime=zeros(length(RatioList),1);
    absDelList=zeros(length(RatioList),2);
    MaxCumsumList=zeros(length(RatioList),2);
    for i=1:length(RatioList)
        tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/tpel.txt']);
    
        temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
        if isempty(temp)
            ejectNo=size(tpel,1);
        else
            ejectNo=temp;
        end
        clear temp;
        ejecttime=tpel(ejectNo,1);
    
        tpel=tpel(1:ejectNo-1,:);
        %absDelList(i,1)=abs(tpel(end,4)-tpel(1,4));
        absDelList(i,1)=max(abs(tpel(:,4)-mean(tpel(:,4))));
        %absDelList(i,2)=abs(tpel(end,3)-tpel(1,3));
        absDelList(i,2)=max(abs(tpel(:,3)-mean(tpel(:,3))));
    
        CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
        ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    
%         di_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',di_name,'.txt']);
        de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);

%         di_record_inout=di_record_inout(1:ejectCENo);
        de_record_inout=de_record_inout(1:ejectCENo);
%         di_sum=cumsum(di_record_inout);
        de_sum=cumsum(de_record_inout);
    
        %EjectTime(i)=ejecttime/365.25;
    
%         MaxCumsumList(i,1)=max(di_sum)-min(di_sum);   
        MaxCumsumList(i,2)=max(de_sum)-min(de_sum);
    end

%     %%% Cum Sum Max
%     [PCumsum1,HCumsum1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)),1);  
%     RCumsum1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),1)));
%     Cumsumfit1=exp(polyval(PCumsum1,log(RatioList1)));
% 
%     [PCumsum2,HCumsum2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)),1);  
%     RCumsum2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,1)));
%     Cumsumfit2=exp(polyval(PCumsum2,log(RatioList2)));

    %%de

    %% Cum Sum Max
    [PCumsumde1,HCumsumde1]=polyfit(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)),1);  
    RCumsumde1=corrcoef(log(RatioList1),log(MaxCumsumList(1:length(RatioList1),2)));
    Cumsumfitde1=exp(polyval(PCumsumde1,log(RatioList1)));

    [PCumsumde2,HCumsumde2]=polyfit(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)),1);  
    RCumsumde2=corrcoef(log(RatioList2),log(MaxCumsumList(length(RatioList1)+1:end,2)));
    Cumsumfitde2=exp(polyval(PCumsumde2,log(RatioList2)));

    %ReviseNo=length(RatioList1)+1;
    
    %%%%%% Cumsum
    subplot(row,col,ipl);
    loglog(RatioList,MaxCumsumList(:,2),'k.');hold all;
    loglog(RatioList,absDelList(:,2),'b.')

    loglog(RatioList1,Cumsumfitde1,'r-','linewidth',linewidth);
    loglog(RatioList2,Cumsumfitde2,'r-','linewidth',linewidth);
    Cumsumfitde1_extend=exp(polyval(PCumsumde1,log(RatioList2)));
    loglog(RatioList2,Cumsumfitde1_extend,'r-.','linewidth',linewidth);

    ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
    xlabel('$m_P/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');

    %if strcmp(PluName,'1999CE119')   
    %    yylim=[1e-7 1e1];
    %else
        yylim=[1e-8 1e0];
    %end
    ylim(yylim);
    xlim(xxlim);

    plot([barrier barrier],yylim,'r--');
    hold off;
    set(gca,'xTick',xxtick);
    
    title(titlename,'fontsize',fontsize);

    set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))));
    set(text(xp,2/10*yylim(2),['$$\ln{y} = ',num2str(PCumsumde1(1),'%.4f'),'\,\ln{x}',num2str(PCumsumde1(2),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize-3,'color','red');
    set(text(xp,0.4/10*yylim(2),['$$R^2 = ',num2str(RCumsumde1(1,2)*RCumsumde1(2,1),'%.4f'),'$$']),...,
    'Interpreter','latex','fontsize',fontsize-3,'color','red');
    %set(gca,'xticklabel',[]);
% set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
% annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
%            '(B2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');

end

%end