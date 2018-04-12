%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname2='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

diname1='di_record_inout';
diname2='di_fit_perturb';
% diname2='di_record_inout';

fontsize=15;
xxlim=[1e-10 10];
yylim=[1e-5 1];
N=100;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    
    subplot(2,2,isub);
    loglog(0,0,'w');hold all;
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
        NCE=length(di_record_inout);
        disp(NCE)

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
    
    patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],...,
        [yylim(1) 1/NCE 1/NCE yylim(1)],...,
        'k','facealpha',0.50,'edgecolor','none');
    
    switch iplot
        case 1
            color='k';
        case 2
            color='r';
    end
    
    Meandi=mean(di_norm);
    plot([Meandi Meandi],yylim,[color,'--'],'linewidth',1.5);

    h=loglog(plotx,ploty,[color,'.-'],'linewidth',1.3,'markersize',10);
    eval(['h',num2str(iplot),'=h;']);
    %area(plotx,ploty,'FaceColor',[0.9 0.9 0.9],'edgecolor','none');

    xlim(xxlim);
    ylim(yylim);

    if iplot==2
        set(text(exp(log(xxlim(1))+(log(xxlim(2))-log(xxlim(1)))/2.5),exp((log(1/NCE)+log(yylim(1)))/2),'\textbf{$$N_{CE}<1$$}'),...,
        'Interpreter','latex','fontsize',fontsize,'color','w');
    end
    
    switch iplot
        case 1; pos=0.2;
        case 2; pos=0.07;
    end
    set(text(5e-10,pos*yylim(2),['$$N= ',num2str(length(di_norm)),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color',color);
    
    end
    legend([h1 h2],{'Test','Mass'},'box','on','Interpreter','latex','fontsize',fontsize/5*4);

    hold off;
end
