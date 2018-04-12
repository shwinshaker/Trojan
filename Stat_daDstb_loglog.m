%%bar di
clear;

ffname='RealPlutinosNpl';
% fname={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};
ffname1='RealPlutinos';
fname={'1999CE119_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172_1Gyr'};

titlename={'1999CE119&2004UP10';'1999CE119&2006RJ103';'2001FU172&2004UP10'};%;'2001FU172&2006RJ103'};

diname='da_record_inout';
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
    

    di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'_40pl_copy/',diname,'.txt']);
    NCE=length(di_record_inout);

    di_norm=di_record_inout;    
    di_norm_0=di_norm(di_norm~=0);
    
    H=kstest2(di_norm_0(di_norm_0>0),-di_norm_0(di_norm_0<0));
    disp(titlename{isub});
    disp(H);
    
    if isub==1
        Minlog=-log(min(abs(di_norm_0)));

        dix=0:N;
        dix=dix';
        dlog=Minlog/N;
        dix=exp(-dix*dlog);
%         dix=[-dix;dix];
        dix=sort(dix);
    end
        
    %%%%%%% Theo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/plel.txt']);
    tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/tpel.txt']);
    r2hill=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/r2hill_record.txt']);

    pla=plel(:,2);plamean=mean(pla);
    ple=plel(:,3);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
    plinc=plel(:,4);plincmean=mean(plinc);%plincstd=std(plinc);
    
    tpa=tpel(:,2);tpamean=mean(tpa);
    tpe=tpel(:,3);tpemean=mean(tpe);
    tpinc=tpel(:,4);tpincmean=mean(tpinc);
    
    N=1000000;%
    %N=length(di_record_inout); %% sample size
    %std0=abs(plincmean)^(1/3)+1;
    switch diname
    case 'di_record_inout'
        [dinc,~,~]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
    case 'de_record_inout'
        [~,dinc,~]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
    case 'da_record_inout'
        [~,~,dinc]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
    end
    dincAbs=abs(dinc);
    
    H=kstest2(di_norm,dinc);
    disp('Theo-Numer test');
    disp(H);
    countx_theo=histcounts(dincAbs,dix)'/length(dincAbs)/2;
    ploty_theo=countx_theo;

    %%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%5%%%%%
    
    countx_p=histcounts(di_norm(di_norm>0),dix)'/NCE;
    countx_n=histcounts(-di_norm(di_norm<0),dix)'/NCE;

    dixx=(dix(2:end)+dix(1:end-1))/2; 
    dixlen=dix(2:end)-dix(1:end-1);

    plotx_p=dixx(dixx>=0);
    plotx=dixx;
    ploty_p=countx_p;
    ploty_n=countx_n;

    subplot(2,2,isub);
    
    loglog(1,1,'w');hold all;
    
    patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],...,
        [yylim(1) 1/NCE 1/NCE yylim(1)],...,
        'k','facealpha',0.75,'edgecolor','none');
    
    absMeandi=abs(mean(di_norm));
    Vardi=var(di_norm);
    plot([absMeandi absMeandi],yylim,'k--');

    Meanabsdi=mean(abs(di_norm));
    Varabsdi=var(abs(di_norm));
    plot([Meanabsdi Meanabsdi],yylim,'k:','linewidth',1.0);

    Sumdi=abs(sum(di_norm));
    plot([Sumdi Sumdi],yylim,'k-.','linewidth',1.0);

    pp=loglog(plotx,ploty_p,'r.','linewidth',1.3,'markersize',10);
    pn=loglog(plotx,ploty_n,'b.','linewidth',1.3,'markersize',10);
    loglog(plotx,ploty_theo,'c-','linewidth',1.3,'markersize',10);

    xlim(xxlim);
    ylim(yylim);

    hold off;
    
    set(text(exp(log(xxlim(1))+(log(xxlim(2))-log(xxlim(1)))/2.5),exp((log(1/NCE)+log(yylim(1)))/2),'\textbf{$$N_{CE}<1$$}'),...,
        'Interpreter','latex','fontsize',fontsize,'color','w');
    
    set(text(xxlim(1)*5,0.3,['$$N_{\Delta ',tag,'>0}= ',num2str(length(find(di_record_inout>=0))),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','r');
    set(text(xxlim(1)*5,0.1,['$$N_{\Delta ',tag,'<0}= ',num2str(length(find(di_record_inout<0))),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','b');
    title(titlename{isub},'fontsize',fontsize);

    set(gca,'xTick',power(10,log10(xxlim(1)):1:log10(xxlim(2))));
    
    xlabel(['$\Delta ',tag,'~',Unit,'$'],'fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end
