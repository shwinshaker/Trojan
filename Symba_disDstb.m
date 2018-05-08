%%bar eh
clear;
ffname='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

Dir='ServerMount';

%% remove points too close
% dishillThresh=0.03;
% di_record_perturb(dishill<dishillThresh)=[];
% dishill(dishill<dishillThresh)=[];

%% scale
hillmax=3.5;
fontsize=15;

% yylim=0.05;
yylim=2.5;
N=50;
dhill=1/N;

figure(1);
set(gcf,'Position',[400,200,650,650/7*5],'color','w');

for isub=1:4
    subplot(2,2,isub);
    fname=fname1{isub};
    titlename=titlename1{isub};
    disp(fname);

    plot(0,0,'w');hold all;
    ierror=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt']);
    PV_record_pl=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt']);
    PV_record_tp=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt']);
    PV_record_pl(ierror,:)=[];
    PV_record_tp(ierror,:)=[];

    CEdis=((PV_record_tp(:,1)-PV_record_pl(:,1)).^2+...,
        (PV_record_tp(:,2)-PV_record_pl(:,2)).^2+(PV_record_tp(:,3)-PV_record_pl(:,3)).^2).^(1/2);
    
    r2hill_record=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt']);
    r2hill_mean=mean(r2hill_record);

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

    [P,H]=polyfit(plotx,ploty,1);
    R=corrcoef(plotx,ploty);
    yfit=polyval(P,plotx);

    bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');hold all;
    plot(plotx,ploty,'k.-','linewidth',1.3,'markersize',10);
    
    plot(plotx,yfit,'r-');
    
    hold off;

    xlim([0 1]);
    ylim([0 yylim]);
    set(text(0.1,8/9*yylim,['$$N_{CE}= ',num2str(length(dishill)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    if P(2)>=0
        set(text(0.1,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x+',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    else
        set(text(0.1,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    end
    set(text(0.1,3.5/5*yylim,['$$R^2=',num2str(R(1,2)*R(2,1),'%.4f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    
%     set(text(0.1,3.0/5*yylim,['$N_{bin}=',num2str(N),'$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    
    title(titlename,'fontsize',fontsize);
%     xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
    xlabel('$\gamma_R$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$PDF$','fontsize',fontsize,'Interpreter','latex');
    
end
