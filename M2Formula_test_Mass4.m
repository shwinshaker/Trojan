

if exist('M2Fit','var') && exist('M2Actual','var')
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')
    clear;
    
%% Normalized fitting -----------------------------------------------
PI=[-0.278135582863818,-0.0344962140121138];
Pe=[-20954.0722195194,131689.157258553,-343200.239196507,474820.375952506,-367870.942458582,151356.608605849,-25840.8838186323];
PTroI=[-0.392834631735606,0.0655236219957504];
Pm=[2.12729693837174,-0.118787510389832];

%% Zero set -----------------------------------------------
de_name='de_record_inout';
Dir='ServerMount';
DDir='MassTest2';
PluName='1999CE119';
fname=[PluName,'_',num2str(sprintf('%.4f',1.0)),'MP'];
tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/tpel.txt']);
plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/plel.txt']);
de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname,'/',de_name,'.txt']);

ple0=mean(plel(:,3));
plI0=mean(plel(:,4));
tpTroI0=mean(tpel(:,4));
Ratio0=1.0;
R0=sum((de_record_inout-mean(de_record_inout)).^2);

fitI=@(I)exp(polyval(PI,sind(I)/sind(plI0)));
fitITro=@(ITro)exp(polyval(PTroI,sind(ITro)/sind(tpTroI0)));
fite=@(e)exp(polyval(Pe,e/ple0));
fitm=@(mp)exp(polyval(Pm,log(mp/Ratio0))); %(mpluto)

%% M2 formula -----------------------------------------------
% M2=@(I,e,mp,ITro)fitI(I)/fitI(I0).*fite(e)/fite(e0).*fitm(mp).*fitI(ITro)/fitI(ITro0);
M2=@(I,e,mp,ITro)R0*fitI(I).*fite(e).*fitm(mp).*fitI(ITro);

%% -----------------------------------
%% Remove Time factor: 4.5 M~NCE*Sigma
%% -----------------------------------

%% G formula
G=@(M2)2*(2/pi*M2).^(1/2);

%% simple test -----------------------------------------------

Gdata=cell(4,2);
M2data=cell(4,2);
Ratiodata=cell(4,1);

Dir='ServerMount';
de_name='de_record_inout';

for ipl=1:4
    
    switch ipl
        case 1
            PluName='1999CE119';
            DDir='MassTest2';
        case 2
            PluName='2001FU172';
            DDir='MassTest3';
        case 3
            PluName='1999CE119_2006RJ103';
            DDir='MassTest4';
        case 4
            PluName='2001FU172_2006RJ103';
            DDir='MassTest5';
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

if strcmp(PluName,'1999CE119') 
    barrier=5.0;
elseif strcmp(PluName,'1999CE119_2006RJ103') 
    barrier=10.0;
else
    barrier=30.0;
end

iF=find(RatioList2<=barrier,1,'last');
RatioList1=[RatioList1;RatioList2(1:iF)];
RatioList2(1:iF)=[];

M2Actual=zeros(length(RatioList1),1);
M2Fit=zeros(length(RatioList1),1);
GActual=zeros(length(RatioList1),1);
GFit=zeros(length(RatioList1),1);

for i=1:length(RatioList1)
    disp(RatioList1(i));
    
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
    
    CE_record=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/CE_record.txt']);
    ejectCENo=find(CE_record(:,1)<ejecttime,1,'last');
    
    de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fname{i},'/',de_name,'.txt']);
    de_record_inout=de_record_inout(1:ejectCENo);  
    
    M2Actual(i)=sum((de_record_inout-mean(de_record_inout)).^2);
    M2Fit(i)=M2(mean(plel(:,4)),mean(plel(:,3)),RatioList1(i),mean(tpel(:,4)));
    
    GFit(i)=G(M2Fit(i));
    de_sum=cumsum(de_record_inout);
    GActual(i)=max(de_sum)-min(de_sum);
end

    Gdata{ipl,1}=GFit;
    Gdata{ipl,2}=GActual;
    M2data{ipl,1}=M2Fit;
    M2data{ipl,2}=M2Actual;
    Ratiodata{ipl}=RatioList1;
    
end

end

figure;
set(gcf,'Position',[400,100,700,500],'color','w');

fontsize=12;

subplot(2,1,1);
loglog(1,1,'w');hold all;
h1=plot(1,1,'r-');
h2=plot(1,1,'m-');
h3=plot(1,1,'b-');
h4=plot(1,1,'c-');
legend([h1 h2 h3 h4],{'1999CE119&2004UP10','2001FU172&2004UP10','1999CE119&2006RJ103','2001FU172&2006RJ103'},...,
    'fontsize',fontsize,'location','best');

for ipl=1:4
    switch ipl
        case 1
            color='r';
        case 2
            color='m';
        case 3
            color='b';
        case 4
            color='c';
    end
    loglog(Ratiodata{ipl},M2data{ipl,2},[color,'.-']);
    loglog(Ratiodata{ipl},M2data{ipl,1},[color,'.:']);
end
xlabel('$m/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$M_{2,e}$','fontsize',fontsize,'Interpreter','latex');

subplot(2,1,2);
loglog(1,1,'w');hold all;
for ipl=1:4
    switch ipl
        case 1
            color='r';
        case 2
            color='m';
        case 3
            color='b';
        case 4
            color='c';
    end
    loglog(Ratiodata{ipl},Gdata{ipl,2},[color,'.-']);
    loglog(Ratiodata{ipl},Gdata{ipl,1},[color,'.:']);
end

xlabel('$m/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
