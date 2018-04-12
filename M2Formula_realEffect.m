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
M2=@(I,e,mp,ITro)fitI(I)/fitI(I0).*fite(e)/fite(e0).*fitm(mp).*fitI(ITro)/fitI(ITro0)*4.5; 
%% Time factor: 4.5 M~NCE*Sigma
G=@(M2)2*(2/pi*M2).^(1/2);

%% test -----------------------------------------------
Dir='ServerMount';
DDir='MassTest3';
fName='2001FU172_0.1496MP';
de_name='de_record_inout';

de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/',de_name,'.txt']);
tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/tpel.txt']);
plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/plel.txt']);

I=plel(1,4);
e=plel(1,3);
mp=0.1496;
TroI=tpel(1,4);

%M2=M20*fitI(I)/fitI(I0)*fite(e)/fite(e0)*fitm(mp)/fitm(m0);
M2Theo=M2(I,e,mp,TroI);
M2Actual=sum((de_record_inout-mean(de_record_inout)).^2);
disp(M2Theo);
disp(M2Actual);

GTheo=G(M2Theo);
de_sum=cumsum(de_record_inout);
GActual=max(de_sum)-min(de_sum);
disp(GTheo);
disp(GActual);

%% test -----------------------------------------------

%% Large Plutinos -----------------------------------------------
largePlutinoName={'Pluto';'Orcus';'2003AZ84';'Ixion';...,
    '2003VS2';'2003UZ413';'2004PF115';'2004UX10';...,
    'Huya';'2006HJ123';'2002XV93';'2001QF298';'1999TC36';...,
    '2002VU130';'2002VR128';'2002VE95'};
MassBase=1e20; %kg
largePlutinoMass=[130;6.32;3;3;1.5;2;3.5;0.3;0.5;0.012;1.7;0.7;0.1275;0.16;1;0.15];
largePlutinoD=[2322;917;727;617;523;600;406.3;361.2;406;283.1;549.2;408.2;393.1;252.9;448.5;249.8];
largePlutinoInc=[17.1;20.6;13.6;19.6;14.8;12.0;13.4;9.5;15.5;12.0;13.3;22.4;8.4;14.0;14.0;16.3];
largePlutinoa=[39.3;39.2;39.4;39.7;39.3;39.2;39.0;39.2;39.4;39.3;39.3;39.3;39.3;39.3;39.3;39.4];
largePlutinoq=[29.7;30.3;32.3;30.1;36.4;30.4;36.5;37.4;28.5;27.4;34.5;34.9;30.6;31.2;28.9;30.4];
largePlutinoe=1-largePlutinoq./largePlutinoa;
largePlutinoMassRatio=largePlutinoMass/largePlutinoMass(1);

%% Total M2
NplTotal=1e6;
MplTotal=1; % 1 pluto mass
Mpl=MplTotal/NplTotal;

Ipl=I0; %% default 1999CE119
epl=e0; %% default 1999CE119

Nbar=100;
MinPower=-6;
MTrojanList=0:Nbar;
MTrojanList=MTrojanList';
MTrojanList=MinPower/100*MTrojanList;
MTrojanList=power(10,MTrojanList);

GTotal=zeros(length(MTrojanList),3);
GLargeMax=zeros(length(MTrojanList),1);
GLeftMax=zeros(length(MTrojanList),1);
GLargeMin=zeros(length(MTrojanList),1);
GLeftMin=zeros(length(MTrojanList),1);

M2Total=zeros(length(MTrojanList),1);
M2LargeTotal=zeros(length(MTrojanList),1);
M2LeftTotal=zeros(length(MTrojanList),1);

ITroMax=0.0;
ITroMin=30.0;
for i=1:length(MTrojanList)    
    MTrojan=MTrojanList(i);
    M2LargeTotal(i,1)=sum(M2(largePlutinoInc,largePlutinoe,largePlutinoMassRatio,ITro0).*((largePlutinoMassRatio+MTrojan)./largePlutinoMassRatio).^2);
    M2LeftTotal(i,1)=NplTotal*M2(Ipl,epl,Mpl,ITro0)*((Mpl+MTrojan)/Mpl)^2;
    M2Total(i,1)=M2LargeTotal(i,1)+M2LeftTotal(i,1);
    GTotal(i,1)=G(M2Total(i,1));

    M2LargeTotalMax=sum(M2(largePlutinoInc,largePlutinoe,largePlutinoMassRatio,ITroMax).*((largePlutinoMassRatio+MTrojan)./largePlutinoMassRatio).^2);
    M2LeftTotalMax=NplTotal*M2(Ipl,epl,Mpl,ITroMax)*((Mpl+MTrojan)/Mpl)^2;
    M2TotalMax=M2LargeTotalMax+M2LeftTotalMax;
    GTotal(i,2)=G(M2TotalMax);
    GLargeMax(i,1)=G(M2LargeTotalMax);
    GLeftMax(i,1)=G(M2LeftTotalMax);
    
    M2LargeTotalMin=sum(M2(largePlutinoInc,largePlutinoe,largePlutinoMassRatio,ITroMin).*((largePlutinoMassRatio+MTrojan)./largePlutinoMassRatio).^2);
    M2LeftTotalMin=NplTotal*M2(Ipl,epl,Mpl,ITroMin)*((Mpl+MTrojan)/Mpl)^2;
    M2TotalMin=M2LargeTotalMin+M2LeftTotalMin;
    GTotal(i,3)=G(M2TotalMin);
    GLargeMin(i,1)=G(M2LargeTotalMin);
    GLeftMin(i,1)=G(M2LeftTotalMin);
end

%% Figure ---------------------------------------------------------------
fontsize=15;
figure;
set(gcf,'Position',[400,100,500,500],'color','w');

loglog(MTrojanList,GTotal(:,1),'w.');hold all;
hr=loglog(MTrojanList,GTotal(:,2),'r-');
hb=loglog(MTrojanList,GTotal(:,3),'b-');

xxlim=get(gca,'xlim');
plot(xxlim,[0.01 0.01],'k--');

indexMax=find(GTotal(:,2)>0.01,1,'last');
indexMin=find(GTotal(:,3)>0.01,1,'last');
MTrojanMax=(MTrojanList(indexMax)+MTrojanList(indexMax+1))/2;
MTrojanMin=(MTrojanList(indexMin)+MTrojanList(indexMin+1))/2;
%yylim=get(gca,'ylim');
yylim=[1e-4 1];
ylim(yylim);
patch([MTrojanMax MTrojanMax MTrojanMin MTrojanMin],[yylim(1) yylim(2) yylim(2) yylim(1)],'y','facealpha',0.3,'edgecolor','none');
patch([xxlim(1) xxlim(1) MTrojanMax MTrojanMax],[yylim(1) yylim(2) yylim(2) yylim(1)],'g','facealpha',0.3,'edgecolor','none');
patch([MTrojanMin MTrojanMin xxlim(2) xxlim(2)],[yylim(1) yylim(2) yylim(2) yylim(1)],'r','facealpha',0.3,'edgecolor','none');

loglog(MTrojanList,GTotal(:,2),'r-');
loglog(MTrojanList,GTotal(:,3),'b-');

loglog(MTrojanList,GLargeMax,'r--');
loglog(MTrojanList,GLeftMax,'r-.');
loglog(MTrojanList,GLargeMin,'b--');
loglog(MTrojanList,GLeftMin,'b-.');

hLarge=loglog(MTrojanList,GLargeMin,'k--');
hLeft=loglog(MTrojanList,GLeftMin,'k-.');

xlabel('$m_T/m_{Pluto}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
legend([hr hb],{['$I_T=',sprintf('%2.1f',ITroMax),'\,\,\rm{DEG}$'],['$I_T=',sprintf('%2.1f',ITroMin),'\,\rm{DEG}$']},'fontsize',fontsize,'location','northwest','Interpreter','latex');

ah=axes('position',get(gca,'position'),'visible','off');
[legh,objh]=legend(ah,[hLarge hLeft],{'LP','ESP'},'fontsize',fontsize/10*8,'location','northwest');
set(legh,'position',[0.167,0.73,0.1,0.08]);

delete([hLarge hLeft]);

hold off;
