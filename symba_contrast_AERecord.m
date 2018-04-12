%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname2='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

plaename='AE_record_pl';
tpaename='AE_record_tp';

fontsize=15;
yylim=0.15;
N=50;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    
    subplot(2,2,isub);
    plot(0,0,'w');hold all;
    % set(gca,'xTick',power(10,-10:1:2));
    title(titlename{isub},'fontsize',fontsize);

    % xlabel('$\Delta I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    % ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');
        
    for iplot=1:2
      
        fname=eval(['fname',num2str(iplot)]);
        ffname=eval(['ffname',num2str(iplot)]);
        
        disp(ffname);
        disp(fname);
        
        AE_record_pl=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',plaename,'.txt']);
        AE_record_tp=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',tpaename,'.txt']);
    
        switch iplot
            case 1
                color='k';
            case 2
                color='r';
        end
        
        h=plot(AE_record_pl(:,1),AE_record_pl(:,2),[color,'.']);
        eval(['h',num2str(iplot),'=h;']);
        
    end
    %legend([h1 h2],{'Rmvs','Symba'},'box','off','Interpreter','latex','fontsize',fontsize/5*4);
    
    hold off;
end
