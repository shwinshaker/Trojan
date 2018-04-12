% CE times - total inc
if exist('CE_times','var') && exist('Totdi','var')
    tag='repeat';
else
    tag='data';
end
if ~strcmp(tag,'repeat')
clear;
Dirname={'pl';'npl'};%;'ran_nplNK'};
pl_name={'1993RO';'1993SB';'1993SC';'1994JR1';'1994TB';'1995HM5';'1995QY9';'1995QZ9';'1996RR20';'1996SZ4';'1996TP66';'1996TQ66';'1997QJ4';'1998HH151';...,
    '1998HK151';'1998HQ151';'1998UR43';'1998US43';'1998WS31';'1998WU31';'1998WV31';'1998WW24';'1998WZ31';'1999CE119';'1999CM158';'1999RK215';'1999TC36';'1999TR11';...,
    '2000CK105';'2000EB173';'2000FB8';'2000FV53';'2000GE147';'2000GN171';'2000YH2';'2001FL194';'2001FR185';'2001FU172';'2001KD77';'2001KN77';'2001KQ77';...,
    '2001KX76';'2001KY76';'2001QF298';'2001QG298';'2001QH298';'2001RU143';'2001RX143';'2001UO18';'2001VN71';'2001YJ140';'2002CE251';'2002CW224';'2002GE32';...,
    '2002GF32';'2002GL32';'2002GV32';'2002GW31';'2002GY32';'2002VD138';'2002VE95';'2002VR128';'2002VX130';'2002XV93';'2003AZ84';'2003FB128';'2003FF128';...,
    '2003FL127';'2003HA57';'2003HD57';'2003HF57';'2003QB91';'2003QH91';'2003QX111';'2003SO317';'2003SR317';'2003TH58';'2003UV292';'2003VS2';'2003WA191';...,
    '2003WU172';'2004DW';'2004EH96';'2004EJ96';'2004EW95';'2004FU148';'2004FW164';'2005EZ296';'2005EZ300';'2005GA187';'2005GB187';'2005GE187';'2005GF187'};

npl_name={'2001KN77_10pl';'2001KN77_20pl';'2001KN77_40pl';...,
    '1999CE119_5pl';'1999CE119_10pl';'1999CE119_15pl';'1999CE119_20pl';...,
    '1999CE119_25pl';'1999CE119_30pl';'1999CE119_40pl';...,
    '2001FU172_5pl';'2001FU172_10pl';'2001FU172_15pl';'2001FU172_20pl';...,
    '2001FU172_25pl';'2001FU172_30pl';'2001FU172_35pl';'2001FU172_40pl';...,
    '1995QY9_10pl';'1998HK151_10pl';'2000FB8_10pl'};

ran_nplNK_name={'1999CE119_10k','1999CE119_20k','1999CE119_30k','1999CE119_40k',...,
    '1999CE119_50k','1999CE119_60k','1999CE119_70k','1999CE119_80k','1999CE119_90k',...,
    '1999CE119_100k','1999CE119_200k','1999CE119_300k',...,
    '1999CE119_400k','1999CE119_500k','1999CE119_600k','1999CE119_700k','1999CE119_800k',...,
    '1999CE119_900k','1999CE119_1M'};
% ,'1999CE119_200k',...,
%     '1999CE119_300k','1999CE119_400k','1999CE119_500k',...,
%     '2001KN77_100k','1999CE119_1M',...,
% };

Ndata=0;
for ii=1:length(Dirname)
    filename=eval([Dirname{ii},'_name']);
    Ndata=Ndata+length(filename);
end
CE_times=cell(Ndata,1);
Totdi=cell(Ndata,1);
Totde=cell(Ndata,1);
name=cell(Ndata,1);
Dir='ServerMount'; % 'swiftdata'
fontsize=15;

di_acc=[];
de_acc=[];

Ndata=0;
for ii=1:length(Dirname)
    disp(Dirname{ii});
    filename=eval([Dirname{ii},'_name']);
    if strcmp(Dirname{ii},'pl')
        ffname='RealPlutinos';
        de_name='de_record_inout';
        di_name='di_record_inout';
    elseif strcmp(Dirname{ii},'ran_nplNK')  
        ffname='RanPlutinosNK';
        de_name='de_fit_perturb';
        di_name='di_fit_perturb';
    elseif strcmp(Dirname{ii},'npl')
        ffname='RealPlutinosNpl';
        de_name='de_record_inout';
        di_name='di_record_inout';
    end
    for i=1:length(filename)
        name{Ndata+i}=[filename{i},'_',Dirname{ii}];
        if strcmp(Dirname{ii},'ran_nplNK')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'npl')
            temp=regexp(char(filename{i}),'_','split');
            fname=[temp{1},'_1Gyr_',temp{2}];
        elseif strcmp(Dirname{ii},'pl')
            fname=[filename{i},'_1Gyr'];
        end
        disp(filename{i});
        di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_name,'.txt'));
        de_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',de_name,'.txt'));

        if ii==1
            di_acc=[di_acc;di_record_perturb];
            de_acc=[de_acc;de_record_perturb];
        end
        
        CE_times{Ndata+i}=size(di_record_perturb,1);
        Totdi{Ndata+i}=abs(sum(di_record_perturb));
        Totde{Ndata+i}=abs(sum(de_record_perturb));
    end
    Ndata=Ndata+length(filename);
end

%% add acc pl
name=[name;'pl_acc'];
CE_times=[CE_times;size(di_acc,1)];
Totdi=[Totdi;abs(sum(di_acc))];
Totde=[Totde;abs(sum(de_acc))];

end

figure;
set(gcf,'Position',[400,100,700,700*0.618],'color','w');
xliml=10^2;
xlimu=10^5;
maxhill=3.5;

plotx=10.^(log10(xliml):log10(xlimu));

for isub=1:2
    subplot(1,2,isub);
    loglog(1,1,'w');hold all;
    for i=1:length(CE_times)
     s=regexp(name{i},'_','split');
     if length(s)==2
         if strcmp(s{2},'pl')
             marker='.';
             markersize=10;
         elseif strcmp(s{2},'acc')
             marker='.';
             markersize=15;
         end
     elseif length(s)==3
         if strcmp(s{3},'npl')
             marker='x';
             markersize=5;
         end
     elseif length(s)==4
             marker='^';
             markersize=5;
     end
        
     switch isub
         case 1
          loglog(CE_times{i},Totdi{i},['k',marker],'Markersize',markersize);
          ylim([1e-4 1e1]);        
          ylabel('$|\sum\Delta I|~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');

          
         case 2
          loglog(CE_times{i},Totde{i},['k',marker],'Markersize',markersize);
          ylim([1e-5 1e-1]);        
          ylabel('$|\sum\Delta e|$','fontsize',fontsize,'Interpreter','latex');
     end

    end

    h1=plot(1,1,'k.');
    h2=plot(1,1,'kx');
%     h3=plot(1,1,'k^');
%     legend([h1 h2 h3],{'$Integ$','$Clone Integ$','$Ran$'},'Interpreter','latex','fontsize',fontsize,'location','northwest');
    legend([h1 h2],{'$Integ$','$Clone Integ$'},'Interpreter','latex','fontsize',fontsize,'location','northwest');

    hold off;
    xlim([xliml xlimu]);
    xlabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
    set(gca,'xTick',power(10,2:1:7));
    
end