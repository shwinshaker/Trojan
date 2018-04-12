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

ITroMin=0.0;
ITroMax=30.0;
Nbar=100;
ITroList=0:Nbar;
ITroList=ITroList';
ITroList=ITroList/Nbar;
ITroList=ITroMin+(ITroMax-ITroMin)*ITroList;

GTotal=zeros(length(ITroList),1);
GLargeTotal=zeros(length(ITroList),1);
GLeftTotal=zeros(length(ITroList),1);
GPluto=zeros(length(ITroList),1);

M2Total=zeros(length(ITroList),1);
M2LargeTotal=zeros(length(ITroList),1);
M2LeftTotal=zeros(length(ITroList),1);
M2Pluto=zeros(length(ITroList),1);

for i=1:length(ITroList)    
    ITro=ITroList(i);
    M2LargeTotal(i)=sum(M2(largePlutinoInc,largePlutinoe,largePlutinoMassRatio,ITro));
    M2LeftTotal(i)=NplTotal*M2(Ipl,epl,Mpl,ITro);
    M2Total(i)=M2LargeTotal(i)+M2LeftTotal(i);
    M2Pluto(i)=M2(largePlutinoInc(1),largePlutinoe(1),largePlutinoMassRatio(1),ITro);
    
    GLargeTotal(i)=G(M2LargeTotal(i));
    GLeftTotal(i)=G(M2LeftTotal(i));
    GTotal(i)=G(M2Total(i));
    GPluto(i)=G(M2Pluto(i));

end

%% Figure ---------------------------------------------------------------
fontsize=15;
figure;
set(gcf,'Position',[400,100,500,500],'color','w');

semilogy(ITroList,GTotal,'w.');hold all;
xlim([ITroMin ITroMax]);
xxlim=get(gca,'xlim');
% plot(xxlim,[0.01 0.01],'k--');
yylim=[1e-7 0.1];
ylim(yylim);

GThresh=0.01;

patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],[yylim(1) GThresh GThresh yylim(1)],'c','facealpha',0.3,'edgecolor','none');
patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],[GThresh yylim(2) yylim(2) GThresh],'r','facealpha',0.3,'edgecolor','none');

htotal=semilogy(ITroList,GTotal,'k-');
hPluto=semilogy(ITroList,GPluto,'r-');
hLarge=semilogy(ITroList,GLargeTotal,'k--');
hLeft=semilogy(ITroList,GLeftTotal,'k-.');

legend([htotal hPluto hLarge hLeft],{'Total','Pluto','LP','ESP'},'fontsize',fontsize/10*8,'location','northeast');
% 
% indexMax=find(GTotal(:,2)>0.01,1,'last');
% indexMin=find(GTotal(:,3)>0.01,1,'last');
% MTrojanMax=(ITroList(indexMax)+ITroList(indexMax+1))/2;
% MTrojanMin=(ITroList(indexMin)+ITroList(indexMin+1))/2;
% %yylim=get(gca,'ylim');

% patch([xxlim(1) xxlim(1) MTrojanMax MTrojanMax],[yylim(1) yylim(2) yylim(2) yylim(1)],'g','facealpha',0.3,'edgecolor','none');
% patch([MTrojanMin MTrojanMin xxlim(2) xxlim(2)],[yylim(1) yylim(2) yylim(2) yylim(1)],'r','facealpha',0.3,'edgecolor','none');
% 
% loglog(ITroList,GTotal(:,2),'r-');
% loglog(ITroList,GTotal(:,3),'b-');
% 
% loglog(ITroList,GLargeMax,'r--');
% loglog(ITroList,GLeftMax,'r-.');
% loglog(ITroList,GLargeMin,'b--');
% loglog(ITroList,GLeftMin,'b-.');
% 
% hLarge=loglog(ITroList,GLargeMin,'k--');
% hLeft=loglog(ITroList,GLeftMin,'k-.');
% 
xlabel('$I_T~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
% legend([hr hb],{['$I_T=',sprintf('%2.1f',ITroMax),'\,\,\rm{DEG}$'],['$I_T=',sprintf('%2.1f',ITroMin),'\,\rm{DEG}$']},'fontsize',fontsize,'location','northwest','Interpreter','latex');
% 
% ah=axes('position',get(gca,'position'),'visible','off');
% [legh,objh]=legend(ah,[hLarge hLeft],{'LP','ESP'},'fontsize',fontsize/10*8,'location','northwest');
% set(legh,'position',[0.167,0.73,0.1,0.08]);
% 
% delete([hLarge hLeft]);
% 
% hold off;
