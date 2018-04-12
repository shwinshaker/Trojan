%%Inc test

clear;

Dir='ServerMount'; 

Npl=1;
DDir1='IncTest';
DDir2='IncTest2';
PluName1='1999CE119';
PluName2='2001FU172';

di_name='di_record_inout';
de_name='de_record_inout';
%IncList=1./exp(log(1):0.1:log(1000));
IncList=0:1:30;
IncList=IncList';

figure;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
plot(0,0,'w');hold all;

for ipl=1:Npl
    
    DDir=eval(['DDir',num2str(ipl)]);
    PluName=eval(['PluName',num2str(ipl)]);
    
for i=1:length(IncList)
    Inc=IncList(i);
    color=[Inc/30.0 0 1-Inc/30.0];
    fname=[PluName,'_',num2str(sprintf('%.1f',Inc)),'Inc'];
    
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',DDir,'/',fname,'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',DDir,'/',fname,'/PV_record_tp.txt'));
    %AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',DDir,'/',fname,'/AE_record_pl.txt'));
    %AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',DDir,'/',fname,'/AE_record_tp.txt'));
    
    CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',DDir,'/',fname,'/CE_record.txt'));
    Time=CE_record(:,1);
    
    VP=PV_record_pl(:,4:6);
    VT=PV_record_tp(:,4:6);
    vdot=VP(:,1).*VT(:,1)+VP(:,2).*VT(:,2)+VP(:,3).*VT(:,3);
    VPnorm=(VP(:,1).^2+VP(:,2).^2+VP(:,3).^2).^(1/2);
    VTnorm=(VT(:,1).^2+VT(:,2).^2+VT(:,3).^2).^(1/2);
    cosphiv=vdot./(VPnorm.*VTnorm);
    
    plot(Time,acosd(cosphiv),'.','color',color);

end


end

hold off;

    
