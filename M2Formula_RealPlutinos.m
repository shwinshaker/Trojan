

if exist('M2Fit','var') && exist('M2Actual','var')
    tag='repeat';
else
    tag='data';
end

if ~strcmp(tag,'repeat')
    clear;
    
% PI=[-8.12809558498434,-12.7654057155199]; %% mean elements fitting
% PI=[-7.7886,-12.7442];
% Pe=[-43593737.9528132,76679583.5924321,-55930842.8159286,21657443.1912317,-4696208.55974271,540788.357742723,-25853.6147173015];
% Pe=[-6.750371507041402e+06,1.279721235683384e+07,-1.005607841976580e+07,4.192952629347689e+06,-9.784494195879544e+05,1.211440103813769e+05,-6.227831519832173e+03];
% PTroI=[-8.48712064336154,-12.6653858795120];
% PTroI=[-8.2484,-12.7845];
% Pm=[2.1273,-12.8497];

%% Normalized fitting
PI=[-0.278135582863818,-0.0344962140121138];
Pe=[-20954.0722195194,131689.157258553,-343200.239196507,474820.375952506,-367870.942458582,151356.608605849,-25840.8838186323];
PTroI=[-0.392834631735606,0.0655236219957504];
Pm=[2.12729693837174,-0.118787510389832];

%% Zero set
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

%% zero set--------------------------------
%m0=mpluto;
%M20=exp(polyval(Pm,log(m0/mpluto))); %1999CE119
% 
% Dir='ServerMount';
% DDir='MassTest2';
% fName='1999CE119_1.00MP';
% 
% plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/plel.txt']);
% tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/tpel.txt']);
% 
% % I0=mean(plel(:,4));
% % e0=mean(plel(:,3));
% 
% I0=plel(1,4);
% e0=plel(1,3);
% ITro0=tpel(1,4);

%% M2 formula -----------------------------------------------
% M2=@(I,e,mp,ITro)fitI(I)/fitI(I0).*fite(e)/fite(e0).*fitm(mp).*fitI(ITro)/fitI(ITro0);
M2=@(I,e,mp,ITro)R0*fitI(I).*fite(e).*fitm(mp).*fitI(ITro);

%% -----------------------------------
%% Remove Time factor: 4.5 M~NCE*Sigma
%% -----------------------------------

%% G formula
G=@(M2)2*(2/pi*M2).^(1/2);

%% simple test -----------------------------------------------
% Dir='ServerMount';
% DDir='MassTest3';
% fName='2001FU172_0.3679MP';
% de_name='de_record_inout';

% de_record_inout=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/',de_name,'.txt']);
% tpel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/tpel.txt']);
% plel=load(['~/Documents/',Dir,'/LAB/CE_realp/',DDir,'/',fName,'/plel.txt']);
% 
% I=mean(plel(:,4));
% e=mean(plel(:,3));
% mp=str2double(cell2mat(regexp(fName,'(?<=_).*(?=M)','match')));
% TroI=mean(tpel(:,4));
% 
% %M2=M20*fitI(I)/fitI(I0)*fite(e)/fite(e0)*fitm(mp)/fitm(m0);
% M2Theo=M2(I,e,mp,TroI);
% M2Actual=sum((de_record_inout-mean(de_record_inout)).^2);
% disp(M2Theo);
% disp(M2Actual);

% GTheo=G(M2Theo);
% de_sum=cumsum(de_record_inout);
% GActual=max(de_sum)-min(de_sum);
% disp(GTheo);
% disp(GActual);

%% test -----------------------------------------------
DDDir='ServerMount';

DDir='RealPlutinos';
Dirname={'1993RO';'1993SB';'1993SC';'1994JR1';'1994TB';'1995HM5';'1995QY9';'1995QZ9';'1996RR20';'1996SZ4';'1996TP66';'1996TQ66';'1997QJ4';'1998HH151';...,
    '1998HK151';'1998HQ151';'1998UR43';'1998US43';'1998WS31';'1998WU31';'1998WV31';'1998WW24';'1998WZ31';'1999CE119';'1999CM158';'1999RK215';'1999TC36';'1999TR11';...,
    '2000CK105';'2000EB173';'2000FB8';'2000FV53';'2000GE147';'2000GN171';'2000YH2';'2001FL194';'2001FR185';'2001FU172';'2001KD77';'2001KN77';'2001KQ77';...,
    '2001KX76';'2001KY76';'2001QF298';'2001QG298';'2001QH298';'2001RU143';'2001RX143';'2001UO18';'2001VN71';'2001YJ140';'2002CE251';'2002CW224';'2002GE32';...,
    '2002GF32';'2002GL32';'2002GV32';'2002GW31';'2002GY32';'2002VD138';'2002VE95';'2002VR128';'2002VX130';'2002XV93';'2003AZ84';'2003FB128';'2003FF128';...,
    '2003FL127';'2003HA57';'2003HD57';'2003HF57';'2003QB91';'2003QH91';'2003QX111';'2003SO317';'2003SR317';'2003TH58';'2003UV292';'2003VS2';'2003WA191';...,
    '2003WU172';'2004DW';'2004EH96';'2004EJ96';'2004EW95';'2004FU148';'2004FW164';'2005EZ296';'2005EZ300';'2005GA187';'2005GB187';'2005GE187';'2005GF187'};

de_name='de_record_inout';

M2Theo=zeros(length(Dirname),1);
M2Actual=zeros(length(Dirname),1);

IList=zeros(length(Dirname),1);
eList=zeros(length(Dirname),1);

for i=1:length(Dirname)
    
        Dir=[Dirname{i},'_1Gyr'];
        disp(Dirname{i});
        
        tpel=load(['~/Documents/',DDDir,'/LAB/CE_realp/',DDir,'/',Dir,'/tpel.txt']);
        plel=load(['~/Documents/',DDDir,'/LAB/CE_realp/',DDir,'/',Dir,'/plel.txt']);

        de_record_inout=load(['~/Documents/',DDDir,'/LAB/CE_realp/',DDir,'/',Dir,'/',de_name,'.txt']);
     
        M2Actual(i)=sum((de_record_inout-mean(de_record_inout)).^2);
        
        %% Fit
        %IList(i)=plel(1,4);
        IList(i)=mean(plel(:,4));
        %eList(i)=plel(1,3);
        eList(i)=mean(plel(:,3));
        mp=1; %% Default: 1 Mpluto
        %TroI=tpel(1,4); %% Default: 2004UP10
        TroI=mean(tpel(:,4));
        
        M2Theo(i)=M2(IList(i),eList(i),mp,TroI);
 
end

end

figure;
plot(1,1,'w.');hold all;

MaxMaxM2=max(max(M2Actual),max(M2Theo));
MaxSize=50;

% SizeActual=zeros(length(M2Actual),1);
% SizeTheo=zeros(length(M2Theo),1);
% for i=1:length(M2Actual)
%     if M2Actual(i)~=0
%         %SizeActual(i)=log(MaxMaxM2)/log(M2Actual(i))*MaxSize;
%         SizeActual(i)=M2Actual(i)/MaxMaxM2*MaxSize;
%     else
%         SizeActual(i)=10;
%     end
% end
% for i=1:length(M2Theo)
%     if M2Theo(i)~=0
%         %SizeTheo(i)=log(MaxMaxM2)/log(M2Theo(i))*MaxSize;
%         SizeTheo(i)=M2Theo(i)/MaxMaxM2*MaxSize;
%     else
%         SizeTheo(i)=10;
%     end
% end
% 
% scatter(IList,eList,SizeTheo,'r','filled');%,'MarkerEdgeColor','w','LineWidth',2);
% scatter(IList,eList,SizeActual,'k');

for i=1:length(IList)
   
    SizeT=M2Theo(i)/MaxMaxM2*MaxSize;
    if SizeT~=0
        plot(IList(i),eList(i),'ro','MarkerSize',SizeT,'MarkerFaceColor','r');
    else
        plot(IList(i),eList(i),'r.');
    end
    
    SizeA=M2Actual(i)/MaxMaxM2*MaxSize;
    if SizeA~=0
        plot(IList(i),eList(i),'ko','MarkerSize',SizeA);
    else
        plot(IList(i),eList(i),'k.');
    end
    
end

% GTheo=G(M2Theo);
% de_sum=cumsum(de_record_inout);
% GActual=max(de_sum)-min(de_sum);
% disp(GTheo);
% disp(GActual);



