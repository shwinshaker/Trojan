%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname2='RealPlutinosAsTp';
ffname0='RealPlutinos';

fname={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};

titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

diname='de_record_inout';
fontsize=15;

switch diname
    case 'di_record_inout'
        xxlim=[1e-10 10];yylim=[1e-5 1];tag='I';Unit='\mathrm{(DEG)}';
    case 'de_record_inout'
        xxlim=[1e-11 1];yylim=[1e-5 1];tag='e';Unit='';
    case 'da_record_inout'
        xxlim=[1e-9 1];yylim=[1e-5 1];tag='a';Unit='\mathrm{(au)}';
end

N=100;

figure;
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:length(fname)
    
    %%%%%%% Theo
    plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname0,'/',fname{isub},'/plel.txt']);
    tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname0,'/',fname{isub},'/tpel.txt']);
    r2hill=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname0,'/',fname{isub},'/r2hill_record.txt']);
    
    aP=mean(plel(:,2));eP=mean(plel(:,3));IP=mean(plel(:,4));
    aT=mean(tpel(:,2));eT=mean(tpel(:,3));IT=mean(tpel(:,4));

    subplot(2,2,isub);
    loglog(1,1,'w');hold all;
    
    for iplot=1:2
        
        switch iplot
            case 1; color='k';suf='_40pl_copy';
            case 2; color='r';suf='_40pl';
%                 if isub==1     
%                     suf='_2';
% 
%                 end
        end

        ffname=eval(['ffname',num2str(iplot)]);
        di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},suf,'/',diname,'.txt']);
        NCE=length(di_record_inout);

        di_norm=abs(di_record_inout);    
        di_norm_0=di_norm(di_norm~=0);
        
        if isub==1 && iplot==1
            Minlog=-log(min(abs(di_norm_0)));
            
            dix=0:N;
            dix=dix';
            dlog=Minlog/N;
            dix=exp(-dix*dlog);
            %         dix=[-dix;dix];
            dix=sort(dix);   
            dixx=(dix(2:end)+dix(1:end-1))/2; 
            dixlen=dix(2:end)-dix(1:end-1);
        end

        %%%%%%%%Theo   
        N=1000000;%
        switch iplot
            case 1; [dinc,de,da]=Fun_diDstb_theo(N,1,aP,eP,IP,aT,eT,IT);
            case 2; [dinc,de,da]=Fun_diDstb_theo_PlAsTp(N,1,aP,eP,IP,aT,eT,IT);
        end
        switch diname
            case 'di_record_inout'
                dd=abs(dinc);
            case 'de_record_inout'
                dd=abs(de);
            case 'da_record_inout'
                dd=abs(da);
        end

        countx_theo=histcounts(dd,dix)'/length(dd);
        ploty_theo=countx_theo;

        countx=histcounts(di_norm,dix)'/NCE;

        plotx=dixx;
        ploty=countx;
        
        patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],...,
            [yylim(1) 1/NCE 1/NCE yylim(1)],...,
            'k','facealpha',0.75,'edgecolor','none');
         
        absMeandi=abs(mean(di_record_inout));
        Vardi=var(di_norm);
        plot([absMeandi absMeandi],yylim,[color,'--']);
        
        Meanabsdi=mean(abs(di_record_inout));
        Varabsdi=var(abs(di_norm));
        plot([Meanabsdi Meanabsdi],yylim,[color,':'],'linewidth',1.0);
        
        Sumdi=abs(sum(di_record_inout));
        plot([Sumdi Sumdi],yylim,[color,'-.'],'linewidth',1.0);
       
        set(text(exp(log(xxlim(1))+(log(xxlim(2))-log(xxlim(1)))/2.5),exp((log(1/NCE)+log(yylim(1)))/2),'\textbf{$$N_{CE}<1$$}'),...,
            'Interpreter','latex','fontsize',fontsize,'color','w');
        set(text(xxlim(1)*5,0.3/iplot^1.4,['$$N = ',num2str(length(di_norm)),'$$']),...,
            'Interpreter','latex','fontsize',fontsize,'color',color);

        loglog(plotx,ploty,[color,'.'],'markersize',10);
        loglog(plotx,ploty_theo,'c-','linewidth',1.3,'markersize',10);

        xlim(xxlim);
        ylim(yylim);
        
        eval(['di_norm_',num2str(iplot),'=di_norm;']);
        
    end
    
    hold off;
    
    H=kstest2(di_norm_1,di_norm_2);
    disp(titlename{isub});
    disp(H);
    
    title(titlename{isub},'fontsize',fontsize);

    set(gca,'xTick',power(10,log10(xxlim(1)):1:log10(xxlim(2))));
    
    xlabel(['$\Delta ',tag,'~',Unit,'$'],'fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end
