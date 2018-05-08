%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname2='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

diname1='di_record_inout';
diname2='di_fit_perturb';

fontsize=15;
yylim=0.10;
N=100;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    
    subplot(2,2,isub);
    semilogx(0,0,'w');hold all;
    set(gca,'xTick',power(10,-10:1:2));
    title(titlename{isub},'fontsize',fontsize);

    xlabel('$\Delta I~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');
        
    for iplot=1:2
      
        fname=eval(['fname',num2str(iplot)]);
        ffname=eval(['ffname',num2str(iplot)]);
        diname=eval(['diname',num2str(iplot)]);
        
        disp(ffname);
        disp(fname);
        
        di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',diname,'.txt']);

%     if isub==1
%         Maxdi=max(abs(di_record_inout));
%     end
    di_norm=abs(di_record_inout);%/Maxdi;
    di_norm_0=di_norm(di_norm~=0);
        
    if isub==1 && iplot==1
        Minlog=-log(min(di_norm_0));

        dix=0:N;
        dix=dix';
        dlog=Minlog/N;
        dix=exp(-dix*dlog);
        dix=sort(dix);
    end

    countx=histcounts(di_norm,dix)'/length(di_norm);

    dixx=(dix(2:end)+dix(1:end-1))/2;
    
    plotx=dixx;
    ploty=countx;
    
    switch iplot
        case 1
            color='k';
        case 2
            color='r';
    end
    
    h=semilogx(plotx,ploty,[color,'.-'],'linewidth',1.3,'markersize',10);
    eval(['h',num2str(iplot),'=h;']);
    %area(plotx,ploty,'FaceColor',[0.9 0.9 0.9],'edgecolor','none');

    xlim([1e-10 10]);
    ylim([0 yylim]);

    Meandi=mean(di_norm);
    plot([Meandi Meandi],[0 yylim],[color,'--'],'linewidth',1.5);

    set(text(5e-10,(9-iplot/1.2)/9*yylim,['$$N= ',num2str(length(di_norm)),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color',color);
    
    end
    legend([h1 h2],{'Test','Mass'},'box','off','Interpreter','latex','fontsize',fontsize/5*4);

    hold off;
end
