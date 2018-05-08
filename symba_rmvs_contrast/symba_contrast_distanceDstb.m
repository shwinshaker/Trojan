clear;
ffname1='RealPlutinosNpl';
ffname2='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

Dir='ServerMount';

%% scale
hillmax=3.5;
fontsize=15;

yylim=0.05;
N=50;
dhill=hillmax/N;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    subplot(2,2,isub);
    plot(0,0,'w');hold all;
    titlename=titlename1{isub};

    for iplot=1:2
        ffname=eval(['ffname',num2str(iplot)]);
        disp(ffname);
        fname=eval(['fname',num2str(iplot),'{isub}']);

        if iplot==1
            color='k';
        else
            color='r';
        end

        PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
        PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));
        
        CEdis=((PV_record_tp(:,1)-PV_record_pl(:,1)).^2+...,
            (PV_record_tp(:,2)-PV_record_pl(:,2)).^2+(PV_record_tp(:,3)-PV_record_pl(:,3)).^2).^(1/2);
        
        r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
        r2hill_mean=mean(r2hill_record);      
        hill=sqrt(r2hill_mean);
        
        dishill=CEdis/hill;
        
        disbar=0:N;
        disbar=disbar';
        disbar=disbar*dhill;
        plotx=disbar(2:end)-dhill/2;
        
        tot=length(dishill);
        count=histcounts(dishill,disbar)';
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
    legend([h1 h2],{'Rmvs','Symba'},'box','off','Interpreter','latex','fontsize',fontsize/5*4,'location','northwest');
    
    xlim([0 3.5]);
    ylim([0 yylim]);
    %     if P(2)>=0
    %         set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x+',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    %     else
    %         set(text(0.3,4/5*yylim,['$$y = ',num2str(P(1),'%.5f'),'\,x',num2str(P(2),'%.5f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    %     end
    %     set(text(0.3,3.5/5*yylim,['$$R^2=',num2str(R(1,2)*R(2,1),'%.4f'),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename,'fontsize',fontsize);
    xlabel('$Dis.~\mathrm{(R_H)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');
end

