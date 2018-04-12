%%bar di
clear;

Dir1='RealPlutinosNpl';
Dir2='RanPlutinosNpl';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_1Gyr_40pl_ran';'2001FU172_1Gyr_40pl_ran';'1999CE119&2006RJ103_1Gyr_40pl_ran';'2001FU172&2006RJ103_1Gyr_40pl_ran'};

titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

%suf='_ran';

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
        Dir=eval(['Dir',num2str(iplot)]);
        diname=eval(['diname',num2str(iplot)]);
        
        disp(Dir);
        disp(fname);
        
        di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',Dir,'/',fname{isub},'/',diname,'.txt']);
        
        %     if isub==1
        %         Maxdi=max(abs(di_record_inout));
        %     end
        di_norm=abs(di_record_inout);%/Maxdi;
        di_norm_0=di_norm(di_norm~=0);
        
        eval(['di_norm_',num2str(iplot),'=di_norm;']);
        
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
        
        AbsSum=abs(sum(di_record_inout));
        plot([AbsSum AbsSum],[0 yylim],[color,'--'],'linewidth',3.0);

    end
    
    set(text(5e-10,(9-1.5/1.2)/9*yylim,['$$N= ',num2str(length(di_norm)),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','r');

    H=kstest2(di_norm_1,di_norm_2,'Alpha',0.1);
    %disp(titlename{isub});
    disp(H);
    
%     set(text(5e-10,(9-1.8)/9*yylim,['$$H_{KS}= ',num2str(H),'$$']),...,
%         'Interpreter','latex','fontsize',fontsize,'color','r');

    
    legend([h1 h2],{'Integ','Ran'},'box','off','Interpreter','latex','fontsize',fontsize/5*4);%

    hold off;
end
