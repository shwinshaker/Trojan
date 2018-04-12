%%bar di
clear;

ffname='RealPlutinosNpl';
ffname1='RealPlutinos';
fname={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};

% fname={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

suf='40pl';

diname='de_record_inout';
fontsize=15;
if strcmp(diname,'di_record_inout')
    xxlim=[1e-11 1];yylim=[1e-5 1];
else
    xxlim=[1e-11 0.1];yylim=[1e-5 1];
end

N=100;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:length(fname)
    
    di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'_',suf,'/',diname,'.txt']);
    NCE=length(di_record_inout);
    
%% eliminate invalid value pour loglog diagram
    
%     if isub==1
%         Maxdi=max(abs(di_record_inout));
%     end
    di_norm=di_record_inout;%/Maxdi;   
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

    %% Theo
    plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/plel.txt']);
    tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/tpel.txt']);
    % r2hill=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/r2hill_record.txt']);

    pla=plel(:,2);plamean=mean(pla);
    ple=plel(:,3);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
    plinc=plel(:,4);plincmean=mean(plinc);%plincstd=std(plinc);
    
    tpa=tpel(:,2);tpamean=mean(tpa);
    tpe=tpel(:,3);tpemean=mean(tpe);
    tpinc=tpel(:,4);tpincmean=mean(tpinc);
    
    N=1000000;%
    %N=length(di_record_inout); %% sample size
    %std0=abs(plincmean)^(1/3)+1;
    if strcmp(diname,'de_record_inout')
        [~,dinc]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
    else
        [dinc]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
    end        
    dincAbs=abs(dinc);
    
    H=kstest2(di_norm,dinc);
    disp('Theo-Numer test');
    disp(H);
    
    %% compare MI
    disp((sum(di_norm-mean(di_norm)).^2));
    disp((sum(dinc-mean(dinc)).^2));
    
    countx_p=histcounts(di_norm(di_norm>0),dix)'/NCE;%length(di_norm(di_norm>0));
    countx_n=histcounts(-di_norm(di_norm<0),dix)'/NCE;%length(-di_norm(di_norm<0));
    countx_theo=histcounts(dincAbs,dix)'/length(dincAbs)/2;

    dixx=(dix(2:end)+dix(1:end-1))/2; 
    dixlen=dix(2:end)-dix(1:end-1);

    plotx_p=dixx(dixx>=0);
    plotx=dixx;
    ploty_p=countx_p;
    ploty_n=countx_n;
    
    ploty_theo=countx_theo;

    %% Plot
    subplot(2,2,isub);
    
    loglog(1,1,'w');hold all;
    
    patch([xxlim(1) xxlim(1) xxlim(2) xxlim(2)],...,
        [yylim(1) 1/NCE 1/NCE yylim(1)],...,
        'k','facealpha',0.70,'edgecolor','none');
    
    absMeandi=abs(mean(di_norm));
    Vardi=var(di_norm);
    plot([absMeandi absMeandi],yylim,'k--');

    Meanabsdi=mean(abs(di_norm));
    Varabsdi=var(abs(di_norm));
    plot([Meanabsdi Meanabsdi],yylim,'k:','linewidth',1.0);

    % TrimMeanabsdi=trimmean(abs(di_norm),1);
    % plot([TrimMeanabsdi TrimMeanabsdi],[ylim(1) ylim(2)],'k--','linewidth',1.3);

    % Rmsabsdi=rms(abs(di_norm));
    % plot([Rmsabsdi Rmsabsdi],[ylim(1) ylim(2)],'k--','linewidth',3.0);

    Sumdi=abs(sum(di_norm));
    plot([Sumdi Sumdi],yylim,'k-.','linewidth',1.0);

    pp=loglog(plotx,ploty_p,'r.','linewidth',1.3,'markersize',10);
    pn=loglog(plotx,ploty_n,'b.','linewidth',1.3,'markersize',10);
    loglog(plotx,ploty_theo,'c-','linewidth',1.3,'markersize',10);
    
    xlim(xxlim);
    ylim(yylim);

    
    set(text(exp(log(xxlim(1))+(log(xxlim(2))-log(xxlim(1)))/2.5),exp((log(1/NCE)+log(yylim(1)))/2),'\textbf{$$N_{CE}<1$$}'),...,
        'Interpreter','latex','fontsize',fontsize,'color','w');

    hold off;
    if strcmp(diname,'di_record_inout')
        tag='I';
    else
        tag='e';
    end
    set(text(xxlim(1)*5,0.2*yylim(2),['$$N_{\Delta ',tag,'>0}= ',num2str(length(find(di_record_inout>=0))),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','r');
    set(text(xxlim(1)*5,0.07*yylim(2),['$$N_{\Delta ',tag,'<0}= ',num2str(length(find(di_record_inout<0))),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','b');
    title(titlename{isub},'fontsize',fontsize);
    %legend([pp pn],{'$\Delta I>0$','$\Delta I<0$'},...,
    %'fontsize',fontsize,'Interpreter','latex','location','northeast','box','off');


    set(gca,'xTick',power(10,log10(xxlim(1)):1:log10(xxlim(2))));
    
    if strcmp(diname,'di_record_inout')
        xlabel(['$\Delta ',tag,'~\mathrm{(DEG)}$'],'fontsize',fontsize,'Interpreter','latex');
    else
        xlabel(['$\Delta ',tag,'$'],'fontsize',fontsize,'Interpreter','latex');
    end
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end
