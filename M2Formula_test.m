clear;

PI=[-7.7886,-12.7442];
Pe=[-6.750371507041402e+06,1.279721235683384e+07,-1.005607841976580e+07,4.192952629347689e+06,-9.784494195879544e+05,1.211440103813769e+05,-6.227831519832173e+03];
Pm=[2.1273,-12.8497];

PITro=[-8.2484,-12.7845];

fitI=@(I)exp(polyval(PI,sind(I)));
fitITro=@(ITro)exp(polyval(PITro,sind(ITro)));
fite=@(e)exp(polyval(Pe,e));
fitm=@(mp)exp(polyval(Pm,log(mp))); %(mpluto)

%% zero set--------------------------------
%m0=mpluto;
%M20=exp(polyval(Pm,log(m0/mpluto))); %1999CE119

Dir='ServerMount';
DDir='MassTest2';
fName='1999CE119_1.00MP';

plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/plel.txt']);
tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/tpel.txt']);

% I0=mean(plel(:,4));
% e0=mean(plel(:,3));

I0=plel(1,4);
e0=plel(1,3);
ITro0=tpel(1,4);

%% M2 formula -----------------------------------------------
M2=@(I,e,mp,ITro)fitI(I)/fitI(I0).*fite(e)/fite(e0).*fitm(mp).*fitI(ITro)/fitI(ITro0);

%% -----------------------------------
%% Remove Time factor: 4.5 M~NCE*Sigma
%% -----------------------------------

%% G formula
G=@(M2)2*(2/pi*M2).^(1/2);

%% test -----------------------------------------------
Dir='ServerMount';

Npl=4;
DDirList={'MassTest2';'MassTest3';'MassTest4';'MassTest5';};
PluNameList={'1999CE119';'2001FU172';'1999CE119_2006RJ103';'2001FU172_2006RJ103';};

RatioList1=1./exp(log(1):0.1:log(1000))';
RatioList1=sort(RatioList1);
RatioList2=exp(log(1):0.1:log(1000))';
RatioList=[RatioList1;RatioList2];

M2Theo=zeros(length(RatioList),1);
M2Actual=zeros(length(RatioList),1);
de_name='de_record_inout';
fname=cell(length(RatioList),1);

Tottime=1e9; % yr

figure;
fontsize=15;
loglog(1,1,'w.');hold all;

for ipl=1:Npl
    
    DDir=DDirList{ipl};
    PluName=PluNameList{ipl};
    
    disp(PluName);
    
    switch ipl
        case 1
            color='r';
        case 2
            color='b';
        case 3
            color='m';
        case 4
            color='c';
    end
    
for i=1:length(RatioList1)
    Ratio=RatioList1(i);
    fname{i}=[PluName,'_',num2str(sprintf('%.4f',Ratio)),'MP'];
end
for i=1:length(RatioList2)
    Ratio=RatioList2(i);
    fname{length(RatioList1)+i}=[PluName,'_',num2str(sprintf('%.2f',Ratio)),'MP'];
end

for i=1:length(RatioList)
    %disp(RatioList);
    
    tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/tpel.txt']);
    plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/plel.txt']);
    
    temp=find(tpel(:,2)>31.0 | tpel(:,2)<29.0,1,'first');
    if isempty(temp)
        ejectNo=size(tpel,1);
    else
        ejectNo=temp;
    end
    clear temp;
    ejecttime=tpel(ejectNo,1)/365;
    
    reviseFactor=Tottime/ejecttime;
    
    tpel=tpel(1:ejectNo-1,:);
    plel=plel(1:ejectNo-1,:);
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)/365<ejecttime,1,'last');
    
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);

    de_record_inout=de_record_inout(1:ejectCENo);
    %de_sum=cumsum(de_record_inout);
        
    %% Fit
    I=plel(1,4);
    e=plel(1,3);
    mp=RatioList(i);
    TroI=tpel(1,4);
    
    M2Theo(i)=M2(I,e,mp,TroI);
    
    M2Actual(i)=sum((de_record_inout-mean(de_record_inout)).^2)*reviseFactor;

end

loglog(RatioList,M2Theo,[color,'-']);
loglog(RatioList,M2Actual,[color,'.--']);
eval(['h',num2str(ipl),'=plot(1,1,[color,''-'']);']);
end
% GTheo=G(M2Theo);
% de_sum=cumsum(de_record_inout);
% GActual=max(de_sum)-min(de_sum);
% disp(GTheo);
% disp(GActual);

xlabel('$m_P/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_2$','fontsize',fontsize,'Interpreter','latex');

legend([h1,h2,h3,h4],{'1999CE119&2004UP10','2001FU172&2004UP10',...,
    '1999CE119&2006RJ103','2001FU172&2006RJ103'},'fontsize',fontsize,'location','best');

% delete(h1);
% delete(h2);
% delete(h3);
% delete(h4);