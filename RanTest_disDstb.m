%%bar eh
clear;
% ffname='RealPlutinosNpl';
% fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
% titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
ffname1='RealPlutinosNpl';
ffname2='RanPlutinosNpl';
fname1={'1999CE119';'2001FU172';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

suf='40pl';

CE_file='CE_record.txt';
%di_file='di_record_inout.txt';
Dir='ServerMount';

% if strcmp(tag,'ran')
%     ierror=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt'));
%     CE_record(ierror,:)=[];
% end


%% remove points too close
% dishillThresh=0.03;
% di_record_perturb(dishill<dishillThresh)=[];
% dishill(dishill<dishillThresh)=[];

%% scale
hillmax=3.5;
fontsize=15;

yylim=2.5;
N=50;
dhill=1/N;%hillmax

figure(1);
set(gcf,'Position',[400,200,650,650/7*5],'color','w');

for isub=1:4
    subplot(2,2,isub);
    plot(0,0,'w');hold all;

    for iplot=1:2
        ffname=eval(['ffname',num2str(iplot)]);
        if iplot==1
            fname=[fname1{isub},'_1Gyr_',suf];
            titlename=titlename1{isub};
            color='k';
        else
            %fname=[fname1{isub},'_40plTimes']; 
            fname=[fname1{isub},'_1Gyr_',suf,'_ran'];
            color='r';
        end

%     CE_record=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',CE_file));
    %di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_file));
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));

    CEdis=((PV_record_tp(:,1)-PV_record_pl(:,1)).^2+...,
        (PV_record_tp(:,2)-PV_record_pl(:,2)).^2+(PV_record_tp(:,3)-PV_record_pl(:,3)).^2).^(1/2);
    
    r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
    r2hill_mean=mean(r2hill_record);

%     hill=sqrt(r2hill_mean);
%     dishill=CEdis/hill;
    hill=sqrt(r2hill_mean);
    Rth=hill*hillmax;
    dishill=CEdis/Rth;
    
    disbar=0:N;
    disbar=disbar';
    disbar=disbar*dhill;
    plotx=disbar(2:end)-dhill/2;

    tot=length(dishill);
    count=histcounts(dishill,disbar)'/dhill;
    ploty=count/tot;

%     [P,H]=polyfit(plotx,ploty,1);
%     R=corrcoef(plotx,ploty);
%     yfit=polyval(P,plotx);

%     bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');hold all;
    h=plot(plotx,ploty,[color,'.-'],'linewidth',1.3,'markersize',10);
    eval(['h',num2str(iplot),'=h;']);

%     plot(plotx,yfit,'r-');
    
    end
    
    hold off;
    legend([h1 h2],{'Integ','Ran'},'box','on','fontsize',fontsize,'location','southeast');

    xlim([0 1]);
    ylim([0 yylim]);
    set(text(0.1,7/9*yylim,['$$N_{CE}= ',num2str(length(dishill)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     if P(2)>=0
%         set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x+',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     else
%         set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     end
%     set(text(0.3,3.5/5*yylim,['$$R^2=',num2str(R(1,2)*R(2,1),'%.4f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename,'fontsize',fontsize);
    xlabel('$\gamma_R$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$PDF$','fontsize',fontsize,'Interpreter','latex');
end

