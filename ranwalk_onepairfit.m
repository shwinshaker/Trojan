% CE times - total inc

Npl=2;

if exist('datadi','var') && exist('datade','var')
    tag='repeat';
else
    tag='data';
end
if ~strcmp(tag,'repeat')
clear;

Npl=2;

%strplot='sum';
strplot='range';

%errP=[-0.6441 -9.3977];

Dirname={'pl';'npl';'ran_nplNK'};

pl_name1={'1999CE119'};
npl_name1={'1999CE119_5pl','1999CE119_10pl','1999CE119_15pl','1999CE119_20pl',...,
    '1999CE119_25pl','1999CE119_30pl','1999CE119_35pl','1999CE119_36pl',...,
    '1999CE119_37pl','1999CE119_38pl','1999CE119_39pl','1999CE119_40pl'};

ran_nplNK_name1={'1999CE119_1k','1999CE119_2k','1999CE119_3k','1999CE119_4k',...,
    '1999CE119_5k','1999CE119_6k','1999CE119_7k','1999CE119_8k','1999CE119_9k',...,
    '1999CE119_10k','1999CE119_20k','1999CE119_30k','1999CE119_40k',...,
    '1999CE119_50k','1999CE119_60k','1999CE119_70k','1999CE119_80k','1999CE119_90k',...,
    '1999CE119_100k'};
% ,'1999CE119_200k','1999CE119_300k',...,
%     '1999CE119_400k','1999CE119_500k','1999CE119_600k','1999CE119_700k','1999CE119_800k',...,
%     '1999CE119_900k','1999CE119_1M'};
ran_npl_name1={};
%    ,'1999CE119_1M'};

pl_name2={'2001FU172'};
npl_name2={'2001FU172_5pl','2001FU172_10pl','2001FU172_15pl','2001FU172_20pl',...,
    '2001FU172_25pl','2001FU172_30pl','2001FU172_35pl','2001FU172_40pl'};

ran_nplNK_name2={'2001FU172_1k','2001FU172_2k','2001FU172_3k','2001FU172_4k',...,
    '2001FU172_5k','2001FU172_6k','2001FU172_7k','2001FU172_8k','2001FU172_9k',...,
    '2001FU172_10k','2001FU172_20k','2001FU172_30k','2001FU172_40k','2001FU172_50k',...,
    '2001FU172_60k','2001FU172_70k','2001FU172_80k','2001FU172_90k',...,
    '2001FU172_100k'};    
% ,'2001FU172_200k','2001FU172_300k',...,
%     '2001FU172_400k','2001FU172_500k','2001FU172_600k','2001FU172_700k','2001FU172_800k',...,
%     '2001FU172_900k','2001FU172_1M'};
ran_npl_name2={};

pl_name3={'1999CE119&2006RJ103'};
npl_name3={};
npl20_name3={};
npl40_name3={'1999CE119&2006RJ103'};
ran_nplNK_name3={'1999CE119&2006RJ103_10k','1999CE119&2006RJ103_20k',...,
    '1999CE119&2006RJ103_30k','1999CE119&2006RJ103_40k',...,
    '1999CE119&2006RJ103_50k','1999CE119&2006RJ103_60k',...,
    '1999CE119&2006RJ103_70k','1999CE119&2006RJ103_80k',...,
    '1999CE119&2006RJ103_90k'};
ran_npl_name3={};

pl_name4={'2001FU172&2006RJ103'};
npl_name4={};
npl20_name4={};
npl40_name4={'2001FU172&2006RJ103'};
ran_nplNK_name4={};
ran_npl_name4={};


Dir='ServerMount'; % 'swiftdata'
fontsize=15;

for ipl=1:Npl
   
   for ii=1:length(Dirname)
    eval([Dirname{ii},'_name=',Dirname{ii},'_name',num2str(ipl)]);
   end

   Ndata=0;
    for ii=1:length(Dirname)
        filename=eval([Dirname{ii},'_name']);
        Ndata=Ndata+length(filename);
    end
    
    CE_times=cell(Ndata,1);
    Totdi=cell(Ndata,1);
    Maxdi=cell(Ndata,1);
    Meanabsdi=cell(Ndata,1);
    Stdabsdi=cell(Ndata,1);
    Meandi=cell(Ndata,1);
    Stddi=cell(Ndata,1);
    MaxTotdi=cell(Ndata,1);
    name=cell(Ndata,1);

    Totde=cell(Ndata,1);
    Maxde=cell(Ndata,1);
    Meanabsde=cell(Ndata,1);
    Stdabsde=cell(Ndata,1);
    Meande=cell(Ndata,1);
    Stdde=cell(Ndata,1);
    MaxTotde=cell(Ndata,1);
    
    Ndata=0;
    di_acc=[];
    de_acc=[];
for ii=1:length(Dirname)
    disp(Dirname{ii});
    filename=eval([Dirname{ii},'_name']);
    if strcmp(Dirname{ii},'pl')
        ffname='RealPlutinos';
        di_name='di_record_inout';
        de_name='de_record_inout';
    elseif strcmp(Dirname{ii},'ran_npl')
        ffname='RanPlutinos';
        di_name='di_fit_perturb';
        de_name='de_fit_perturb';
    elseif strcmp(Dirname{ii},'ran_nplNK')
        ffname='RanPlutinosNK';
        di_name='di_fit_perturb';
        de_name='de_fit_perturb';
    elseif strcmp(Dirname{ii},'npl')
        ffname='RealPlutinosNpl';
        di_name='di_record_inout';
        de_name='de_record_inout';
    end
    for i=1:length(filename)
        name{Ndata+i}=[filename{i},'_',Dirname{ii}];
        if strcmp(Dirname{ii},'ran_npl')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'ran_nplNK')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'npl')
            temp=regexp(char(filename{i}),'_','split');
            fname=[temp{1},'_1Gyr_',temp{2}];
        else
            fname=[filename{i},'_1Gyr'];
        end
        disp(filename{i});
        di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_name,'.txt'));
        de_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',de_name,'.txt'));

%         p_pick=di_record_perturb>0;
%         n_pick=di_record_perturb<0;
%         di_record_perturb_p=di_record_perturb(p_pick);
%         di_record_perturb_n=di_record_perturb(n_pick);
%         di_record_perturb_p_revise=di_record_perturb_p.*(1+exp(polyval(errP,log(abs(di_record_perturb_p)))));
%         di_record_perturb_n_revise=di_record_perturb_n.*(1-exp(polyval(errP,log(abs(di_record_perturb_n)))));
%         di_record_perturb(p_pick)=di_record_perturb_p_revise;
%         di_record_perturb(n_pick)=di_record_perturb_n_revise;
        %Acc di
        di_acc=[di_acc;di_record_perturb];
        de_acc=[de_acc;de_record_perturb];
        
        if exist(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'),'file')
            r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
            r2hill_mean=mean(r2hill_record);
            hill=sqrt(r2hill_mean);
            disp(hill);
        end

        CE_times{Ndata+i}=size(di_record_perturb,1);
        di_sum=cumsum(di_record_perturb);
        de_sum=cumsum(de_record_perturb);
        
        Totdi{Ndata+i}=abs(sum(di_record_perturb));
        Totde{Ndata+i}=abs(sum(de_record_perturb));

        if CE_times{Ndata+i}==0
         Maxdi{Ndata+i}=0;
         Meanabsdi{Ndata+i}=0;
         Meandi{Ndata+i}=0;
         Stddi{Ndata+i}=0;
         MaxTotdi{Ndata+i}=0;
         Stdabsdi{Ndata+i}=0;
         
         Maxde{Ndata+i}=0;
         Meanabsde{Ndata+i}=0;
         Meande{Ndata+i}=0;
         Stdde{Ndata+i}=0;
         MaxTotde{Ndata+i}=0;
         Stdabsde{Ndata+i}=0;
        else
         Maxdi{Ndata+i}=max(abs(di_record_perturb));
         Meanabsdi{Ndata+i}=mean(abs(di_record_perturb));
         Meandi{Ndata+i}=abs(mean(di_record_perturb));
         Stddi{Ndata+i}=std(di_record_perturb);
         Stdabsdi{Ndata+i}=std(abs(di_record_perturb));
         MaxTotdi{Ndata+i}=max(di_sum)-min(di_sum);
         
         Maxde{Ndata+i}=max(abs(de_record_perturb));
         Meanabsde{Ndata+i}=mean(abs(de_record_perturb));
         Meande{Ndata+i}=abs(mean(de_record_perturb));
         Stdde{Ndata+i}=std(de_record_perturb);
         Stdabsde{Ndata+i}=std(abs(de_record_perturb));
         MaxTotde{Ndata+i}=max(de_sum)-min(de_sum);

         end
    end
    Ndata=Ndata+length(filename);
end

namedatadi=[name CE_times Maxdi Totdi MaxTotdi Meanabsdi Stdabsdi Meandi Stddi];
namedatadi=sortrows(namedatadi,-2);
datadi=cell2mat(namedatadi(:,2:end));
StdFitdi=std(di_acc);
eval(['namedatadi',num2str(ipl),'=namedatadi;']);
eval(['datadi',num2str(ipl),'=datadi;']);
eval(['StdFitdi',num2str(ipl),'=StdFitdi;']);

namedatade=[name CE_times Maxde Totde MaxTotde Meanabsde Stdabsde Meande Stdde];
namedatade=sortrows(namedatade,-2);
datade=cell2mat(namedatade(:,2:end));
StdFitde=std(de_acc);
eval(['namedatade',num2str(ipl),'=namedatade;']);
eval(['datade',num2str(ipl),'=datade;']);
eval(['StdFitde',num2str(ipl),'=StdFitde;']);

end
end

figure;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.4;

row=4;
col=2;
set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
xliml=10^2;
xlimu=10^6;

ylimudi=10; %DEG
ylimldi=1e-4;

ylimude=1;
ylimlde=1e-5;

maxhill=3.5;

plotx=10.^(log10(xliml):0.01:log10(xlimu));

for isub=1:6
    
    for ipl=1:Npl
        switch ipl
            case 1
                color='r';%[255 61 53]/255;
            case 2
                color='b';%[40 160 255]/255;
            case 3
                color='g';
            case 4
                color='c';
        end
        datadi=eval(['datadi',num2str(ipl)]);
        datade=eval(['datade',num2str(ipl)]);

%         [P,H]=polyfit(log(data(:,1)),log(data(:,8)),0);
%          R=corrcoef(log(data(:,1)),log(data(:,8)));
%          StdFit=exp(polyval(P,log(plotx)));
        StdFitdi=eval(['StdFitdi',num2str(ipl)]);
        StdFitde=eval(['StdFitde',num2str(ipl)]);
        
        namedatadi=eval(['namedatadi',num2str(ipl)]);        
        namedatade=eval(['namedatade',num2str(ipl)]);

        subplot(row,col,isub);
        loglog(1,1,'w');hold all;   
        
    for i=1:size(datadi,1)
         s=regexp(char(namedatadi(i,1)),'_','split');
        if length(s)==2
         if strcmp(s{2},'pl')
             marker='o';
             markersize=5;
         elseif strcmp(s{2},'acc')
             marker='o';
             markersize=5;
         end
        elseif length(s)==3
            if strcmp(s{3},'npl')
             marker='s';
             markersize=5;
            end
        elseif length(s)==4
             marker='^';
             markersize=5;
        end
        
     switch isub
         case 1
             loglog(datadi(i,1),datadi(i,3),[color,marker],'Markersize',markersize);
             Efit=(2/pi)^(1/2)*datadi(i,1).^(1/2).*datadi(i,8);
             plot(datadi(i,1),Efit,'marker',marker,'color',color,'Markersize',markersize,'markerfacecolor',color);    
             plot([datadi(i,1) datadi(i,1)],[datadi(i,3) Efit],'-','color',color);
            
         case 3
             loglog(datadi(i,1),datadi(i,4),[color,marker],'Markersize',markersize);
             Efit=2*(2/pi)^(1/2)*datadi(i,1).^(1/2).*datadi(i,8);
             plot(datadi(i,1),Efit,'marker',marker,'color',color,'Markersize',markersize,'markerfacecolor',color);    
             plot([datadi(i,1) datadi(i,1)],[datadi(i,4) Efit],'-','color',color);
             

         case 5
             loglog(datadi(i,1),datadi(i,8),[color,marker],'Markersize',markersize);
             
         case 2
             loglog(datade(i,1),datade(i,3),[color,marker],'Markersize',markersize);
             Efit=(2/pi)^(1/2)*datade(i,1).^(1/2).*datade(i,8);
             plot(datade(i,1),Efit,'marker',marker,'color',color,'Markersize',markersize,'markerfacecolor',color);    
             plot([datade(i,1) datade(i,1)],[datade(i,3) Efit],'-','color',color);
            
         case 4
             loglog(datade(i,1),datade(i,4),[color,marker],'Markersize',markersize);
             Efit=2*(2/pi)^(1/2)*datade(i,1).^(1/2).*datade(i,8);
             plot(datade(i,1),Efit,'marker',marker,'color',color,'Markersize',markersize,'markerfacecolor',color);    
             plot([datade(i,1) datade(i,1)],[datade(i,4) Efit],'-','color',color);
             

         case 6
             loglog(datade(i,1),datade(i,8),[color,marker],'Markersize',markersize);           
     end
            
     end
      
    
    switch isub
        case 1
            loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*StdFitdi,'-','color',color);
            h1=plot(1,1,'r-');
            h2=plot(1,1,'b-');
            legend([h1 h2],{'1999CE119&2004UP10','2001FU172&2004UP10'},'fontsize',fontsize/10*7);

        case 3
            loglog(plotx,2*(2/pi)^(1/2)*plotx.^(1/2).*StdFitdi,'-','color',color);
            
        case 5
            loglog([xliml xlimu],[StdFitdi StdFitdi],'-','color',color);
            
        case 2
            loglog(plotx,(2/pi)^(1/2)*plotx.^(1/2).*StdFitde,'-','color',color);

        case 4
            loglog(plotx,2*(2/pi)^(1/2)*plotx.^(1/2).*StdFitde,'-','color',color);
            
        case 6
            loglog([xliml xlimu],[StdFitde StdFitde],'-','color',color);
    end
        

    end
    
    hold off;
    xlim([xliml xlimu]);
    set(gca,'xTick',power(10,2:1:7));
    switch isub
        case 1
            ylabel('$|A_I|=|\sum\Delta I|~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
            
            set(gca,'xticklabel',[]);
            set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
            annotation('textbox',[Width+LeftRetainWidth/2 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
            
            yylim=[1e-3 1e1]; 
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));

        case 3
            ylabel('$G_I$','fontsize',fontsize,'Interpreter','latex');
             
            set(gca,'xticklabel',[]);
            set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
            annotation('textbox',[Width+LeftRetainWidth/2 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(B1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
       
            yylim=[1e-3 1e1];
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

        case 5
            ylabel('$s_I$','fontsize',fontsize,'Interpreter','latex');
            
            set(gca,'position',[LeftRetainWidth BottomRetainWidth+Height Width Height]);
            annotation('textbox',[Width+LeftRetainWidth/2 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(C1)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
       
            yylim=[1e-5 1e-2];
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1)):1:log10(yylim(2))-1));
            set(gca,'yticklabel',{'10^{-5}','10^{-4}','10^{-3}'});
            xlabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

            
        case 2
            ylabel('$|A_e|=|\sum\Delta e|$','fontsize',fontsize,'Interpreter','latex');
            
            set(gca,'xticklabel',[]);
            set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
            annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+3*Height 0.05 0.05],'edgecolor','none','string',...,
           '(A2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
            
            yylim=[1e-5 1e-1]; 
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))));

        case 4
            ylabel('$G_e$','fontsize',fontsize,'Interpreter','latex');
             
            set(gca,'xticklabel',[]);
            set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
            annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+2*Height 0.05 0.05],'edgecolor','none','string',...,
           '(B2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
       
            yylim=[1e-5 1e-1];
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1))+1:1:log10(yylim(2))-1));

        case 6
            ylabel('$s_e$','fontsize',fontsize,'Interpreter','latex');
            
            set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+Height Width Height]);
            annotation('textbox',[2*Width+LeftRetainWidth*1.5 0.025+Height 0.05 0.05],'edgecolor','none','string',...,
           '(C2)','fontweight','bold','fontsize',fontsize/10*8,'color','k');
       
            yylim=[1e-6 1e-3];
            ylim(yylim);
            set(gca,'yTick',power(10,log10(yylim(1)):1:log10(yylim(2))-1));
            set(gca,'yticklabel',{'10^{-6}','10^{-5}','10^{-4}'});
            xlabel('$N_{CE}$','fontsize',fontsize,'Interpreter','latex');

    end
end