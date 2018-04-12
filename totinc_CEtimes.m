% CE times - total inc
if exist('data','var') && exist('namedata','var')
    tag='repeat';
else
    tag='data';
end
if ~strcmp(tag,'repeat')
clear;
Dirname={'pl';'npl';'npl20';'npl40';'ran_npl'};
%ffname='RealPlutinos';
% filename=['~/Documents/',Dir,'/LAB/tools/Plutino/PlutinoList.txt'];
% fid=fopen(filename,'r');
% Pluinfo=textscan(fid,'%s %f %f %f %f %f %f %f %f %f');
% fclose(fid);
% PluPhi=Pluinfo{10};
% pl_name=unique(Pluinfo{:,1});
pl_name={'1993RO';'1993SB';'1993SC';'1994JR1';'1994TB';'1995HM5';'1995QY9';'1995QZ9';'1996RR20';'1996SZ4';'1996TP66';'1996TQ66';'1997QJ4';'1998HH151';...,
    '1998HK151';'1998HQ151';'1998UR43';'1998US43';'1998WS31';'1998WU31';'1998WV31';'1998WW24';'1998WZ31';'1999CE119';'1999CM158';'1999RK215';'1999TC36';'1999TR11';...,
    '2000CK105';'2000EB173';'2000FB8';'2000FV53';'2000GE147';'2000GN171';'2000YH2';'2001FL194';'2001FR185';'2001FU172';'2001KD77';'2001KN77';'2001KQ77';...,
    '2001KX76';'2001KY76';'2001QF298';'2001QG298';'2001QH298';'2001RU143';'2001RX143';'2001UO18';'2001VN71';'2001YJ140';'2002CE251';'2002CW224';'2002GE32';...,
    '2002GF32';'2002GL32';'2002GV32';'2002GW31';'2002GY32';'2002VD138';'2002VE95';'2002VR128';'2002VX130';'2002XV93';'2003AZ84';'2003FB128';'2003FF128';...,
    '2003FL127';'2003HA57';'2003HD57';'2003HF57';'2003QB91';'2003QH91';'2003QX111';'2003SO317';'2003SR317';'2003TH58';'2003UV292';'2003VS2';'2003WA191';...,
    '2003WU172';'2004DW';'2004EH96';'2004EJ96';'2004EW95';'2004FU148';'2004FW164';'2005EZ296';'2005EZ300';'2005GA187';'2005GB187';'2005GE187';'2005GF187'};
% pl_name={'2001KN77';'1999CE119';'1998WV31';'2001FU172';'1993SB'; ...,
%     '2001UO18';'1994TB';'1999CM158';'2000YH2';'2002CE251';'2000FB8';'1996TP66';'1998HQ151';...,
%     '2001RX143';'1995QY9';'1996SZ4';...,
%     '2000GE147';'2001KD77';'2001FR185';'1998HK151';...,
%     '2002CW224';'2001KY76';'2000CK105';'1995HM5'};

%ffname_npl='RealPlutinosNpl';
npl_name={'2001KN77';'1999CE119';'1995QY9';'1998HK151';'2000FB8'};
npl20_name={'1999CE119','2001KN77'};
npl40_name={'1999CE119','2001KN77','2001FU172','1999CE119&2006RJ103','2001FU172&2006RJ103'};

%ranffname='RanPlutinos';
ran_npl_name={'1999CE119_10k','1999CE119_100k','1999CE119_200k','1999CE119_300k','1999CE119_400k','1999CE119_500k','2001KN77_100k','1999CE119_1M',...,
    '1999CE119_100k_30000Fchk','2004UP10_100k_WideRange'};
%Ndata=length(pl_name)+length(npl_name)+length(npl20_name)+length(npl40_name)+length(ran_npl_name);
Ndata=0;
for ii=1:length(Dirname)
    filename=eval([Dirname{ii},'_name']);
    Ndata=Ndata+length(filename);
end
CE_times=cell(Ndata,1);
Totdi=cell(Ndata,1);
Maxdi=cell(Ndata,1);
Meanabsdi=cell(Ndata,1);
Varabsdi=cell(Ndata,1);
Meandi=cell(Ndata,1);
Vardi=cell(Ndata,1);
MaxTotdi=cell(Ndata,1);
name=cell(Ndata,1);
% npl_name_suffix=cell(length(npl_name),1);
% npl20_name_suffix=cell(length(npl20_name),1);
% npl40_name_suffix=cell(length(npl40_name),1);
% ran_npl_name_suffix=cell(length(ran_npl_name),1);
Dir='ServerMount'; % 'swiftdata'
fontsize=15;

di_acc=[];
Ndata=0;
for ii=1:length(Dirname)
    disp(Dirname{ii});
    filename=eval([Dirname{ii},'_name']);
    if strcmp(Dirname{ii},'pl')
        ffname='RealPlutinos';
        %di_name='de_record_perturb';
        di_name='de_record_inout';
        CE_name='CE_record';
    elseif strcmp(Dirname{ii},'ran_npl')
        ffname='RanPlutinos';
        di_name='de_fit_perturb';
        CE_name='ran_record';
    elseif intersect(Dirname{ii},['npl' 'npl20' 'npl40'])
        ffname='RealPlutinosNpl';
        %di_name='de_record_perturb';
        di_name='de_record_inout';
        CE_name='CE_record';
    end
    for i=1:length(filename)
        name{Ndata+i}=[filename{i},'_',Dirname{ii}];
        if strcmp(Dirname{ii},'ran_npl')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'npl')
            fname=[filename{i},'_1Gyr_10pl'];
            %fname=[filename{i},'_1Gyr'];           
        elseif strcmp(Dirname{ii},'npl20')
            fname=[filename{i},'_1Gyr_20pl'];
        elseif strcmp(Dirname{ii},'npl40')
            fname=[filename{i},'_1Gyr_40pl'];
        else
            fname=[filename{i},'_1Gyr'];
        end
        disp(filename{i});
        di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_name,'.txt'));
        %plel=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/plel.txt'));

        if ii==1
            di_acc=[di_acc;di_record_perturb];
        end
        CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',CE_name,'.txt'));
        if exist(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'),'file')
            r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
            r2hill_mean=mean(r2hill_record);
            hill=sqrt(r2hill_mean);
            disp(hill);
        end


        CE_times{Ndata+i}=size(di_record_perturb,1);
        di_sum=cumsum(di_record_perturb);
        Totdi{Ndata+i}=abs(sum(di_record_perturb));
        if CE_times{Ndata+i}==0
         Maxdi{Ndata+i}=0;
         Meanabsdi{Ndata+i}=0;
         Meandi{Ndata+i}=0;
         Vardi{Ndata+i}=0;
         MaxTotdi{Ndata+i}=0;
         Varabsdi{Ndata+i}=0;
        else
         Maxdi{Ndata+i}=max(abs(di_record_perturb));
         Meanabsdi{Ndata+i}=mean(abs(di_record_perturb));
         Meandi{Ndata+i}=abs(mean(di_record_perturb));
         Vardi{Ndata+i}=var(di_record_perturb);
         Varabsdi{Ndata+i}=var(abs(di_record_perturb));
         MaxTotdi{Ndata+i}=max(di_sum)-min(di_sum);%max(abs(di_sum));
         
         disp(min(CE_record(:,2)));
        end
    end
    Ndata=Ndata+length(filename);
end

%% add acc pl
name=[name;'pl_acc'];
CE_times=[CE_times;size(di_acc,1)];
Maxdi=[Maxdi;max(abs(di_acc))];
Totdi=[Totdi;abs(sum(di_acc))];
di_sum=cumsum(di_acc);
MaxTotdi=[MaxTotdi;max(di_sum)-min(di_sum)];
Meanabsdi=[Meanabsdi;mean(abs(di_acc))];
Varabsdi=[Varabsdi;var(abs(di_acc))];
Meandi=[Meandi;abs(mean(di_acc))];
Vardi=[Vardi;var(di_acc)];

namedata=[name CE_times Maxdi Totdi MaxTotdi Meanabsdi Varabsdi Meandi Vardi];
namedata=sortrows(namedata,-2);
data=cell2mat(namedata(:,2:end));
maxdi=max([max(data(:,2)) max(data(:,3)) max(data(:,4))])+0.02;
end

figure;
set(gcf,'Position',[400,100,1200,500],'color','w');
xliml=10^2;
xlimu=10^7;
%ylimu=10^2; %DEG
%yliml=1e-7;
maxhill=3.5;

plotx=10.^(log10(xliml):log10(xlimu));

for isub=1:4
    %iplot=isub+1;
    subplot(1,4,isub);
    %semilogx(0,0,'w');hold all;
    loglog(1,1,'w');hold all;
    for i=1:length(data)
     %semilogx(data(:,1),data(:,iplot),'k--');
     %loglog(data(:,1),data(:,iplot),'k--');
     s=regexp(char(namedata(i,1)),'_','split');
     if length(s)==2
         if strcmp(s{2},'pl')
             marker='.';
             markersize=10;
         elseif strcmp(s{2},'npl')
             marker='x';
             markersize=5;
         elseif strcmp(s{2},'npl20')
             marker='x';
             markersize=5;
         elseif strcmp(s{2},'npl40')
             marker='x';
             markersize=5;
         elseif strcmp(s{2},'acc')
             marker='.';
             markersize=15;
         end
     elseif length(s)==4
             marker='^';
             markersize=5;
     end
     
     if isub==1
        loglog(data(i,1),data(i,2),['k',marker],'Markersize',markersize);
        %loglog(data(i,1),data(i,7),['r',marker],'Markersize',markersize);
        errorbar(data(i,1),data(i,7),data(i,8),['r',marker],'clipping','off');
        errorbar(data(i,1),data(i,5),2e2*data(i,6),['b',marker],'clipping','off');
        ylim([1e-8 1e2]);
     end
     
     if isub==1 && i==1
        dismin=1./plotx/0.0114;
        
        
        miu=1.94148e-12;
        mius=(2*pi)^2/(365.25)^2;
        at=30.0;
        ap=40.0;
        Rnorm=30.0;
        C=2*miu*Rnorm/(mius*at)^(1/2);
        
        coswfmax=1.0;
        coswfmin=0.3;
        Vt=(mius*(2/Rnorm-1/at))^(1/2);
        Vp=(mius*(2/Rnorm-1/ap))^(1/2);
        Vrnormmax=(Vt^2+Vp^2-2*Vt*Vp*cosd(90))^(1/2);
        Vrnormmin=(Vt^2+Vp^2-2*Vt*Vp*cosd(0))^(1/2);
        Cmax=C/Vrnormmin*coswfmax;
        Cmin=C/Vrnormmax*coswfmin;
        
        fitMaxdiCore=(1-dismin.^2/maxhill^2).^(1/2)*1./(dismin*hill)/pi*180;
        loglog(plotx,Cmax*fitMaxdiCore,'r-');
        loglog(plotx,Cmin*fitMaxdiCore,'r--');
     end
        
     if isub==2
        loglog(data(i,1),data(i,3),['k',marker],'Markersize',markersize);
        Efit=(2/pi)^(1/2)*data(i,1).^(1/2).*data(i,8).^(1/2);
        %errorbar(data(i,1),(2/pi)^(1/2)*data(i,1).^(1/2).*data(i,8).^(1/2),(1-2/pi)*data(i,1).*data(i,8),['r',marker],'clipping','off');
        plot(data(i,1),Efit,['r',marker],'Markersize',markersize);    
        plot([data(i,1) data(i,1)],[data(i,3) Efit],'r-')
        ylim([1e-5 1e2]);        
     end
     if isub==3
        loglog(data(i,1),data(i,4),['k',marker],'Markersize',markersize);         
        Efit=(2/pi)^(1/2)*data(i,1).^(1/2).*data(i,8).^(1/2)*2.0;
        %errorbar(data(i,1),(2/pi)^(1/2)*data(i,1).^(1/2).*data(i,8).^(1/2),(1-2/pi)*data(i,1).*data(i,8),['r',marker],'clipping','off');
        plot(data(i,1),Efit,['r',marker],'Markersize',markersize);    
        plot([data(i,1) data(i,1)],[data(i,4) Efit],'r-')
        ylim([1e-5 1e2]);
     end
     
     if isub==4
        loglog(data(i,1),sqrt(data(i,8)),['k',marker],'Markersize',markersize);  
        ylim([1e-5 1e-2]);
     end
    
     %loglog(data(i,1),data(i,5),['b',marker],'Markersize',markersize);
     %loglog(data(i,1),data(i,6),['r',marker],'Markersize',markersize);
    end
    if isub==2
        %loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*data(8,8).^(1/2),'r--');
        loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*var(di_acc).^(1/2),'r-');
    end
    if isub==3
        %loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*data(8,8).^(1/2)*2.0,'r--');
        loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*var(di_acc).^(1/2)*2.0,'r-');
    end
    hold off;
    xlim([xliml xlimu]);
    xlabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');
    set(gca,'xTick',power(10,2:1:7));
    switch isub
        case 1
            ylabel('$Max~|\Delta Inc|~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
        case 2
            ylabel('$|Sum~\Delta Inc|~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
        case 3
            %ylabel('Max |Accumalative \deltainc| /DEG','fontsize',fontsize);
            ylabel('$Range~\Delta Inc~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');

    end
end

% subplot(1,3,2);
% semilogx(0,0,'w');hold all;
% for i=1:length(data)
%     s=regexp(char(namedata(i,1)),'_','split');
%     if length(s)==1
%         semilogx(data(i,1),data(i,2),'k+');
%     elseif length(s)==2
%         if strcmp(s{2},'10pl')
%             semilogx(data(i,1),data(i,2),'k^');
%         elseif strcmp(s{2},'20pl')
%             semilogx(data(i,1),data(i,2),'k*');
%         elseif strcmp(s{2},'40pl')
%             semilogx(data(i,1),data(i,2),'kx');
%         end
%     elseif length(s)==3
%             semilogx(data(i,1),data(i,2),'rx');
%     end
% end
% xlim([xlim1 10^6]);
% ylim([0 maxdi]);
% xlabel('CE times','fontsize',fontsize);
% ylabel('|Final Accumalative \deltainc| /DEG','fontsize',fontsize);
% 
% subplot(1,3,3);
% semilogx(0,0,'w');hold all;
% for i=1:length(data)
%     s=regexp(char(namedata(i,1)),'_','split');
%     if length(s)==1
%         semilogx(data(i,1),data(i,4),'k+');
%     elseif length(s)==2
%         if strcmp(s{2},'10pl')
%             semilogx(data(i,1),data(i,4),'k^');
%         elseif strcmp(s{2},'20pl')
%             semilogx(data(i,1),data(i,4),'k*');
%         elseif strcmp(s{2},'40pl')
%             semilogx(data(i,1),data(i,4),'kx');
%         end
%    elseif length(s)==3
%             semilogx(data(i,1),data(i,4),'rx');    
%     end
% end
% xlim([xlim1 10^6]);
% ylim([0 maxdi]);
% xlabel('CE times','fontsize',fontsize);
% ylabel('Max |Accumalative \deltainc| /DEG','fontsize',fontsize);