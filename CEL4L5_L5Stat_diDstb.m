%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname={'RealTrojansNpl';'RealPlutinosNpl'};
fname={'2013KY18_1Gyr_40pl';'1999CE119_1Gyr_40pl'};
%fname={'2008LC18_1Gyr_40pl';'1999CE119_1Gyr_40pl'};


titlename={'1999CE119&2013KY18';'1999CE119&2004UP10'};
%titlename={'1999CE119&2008LC18';'1999CE119&2004UP10'};

Nf=length(fname);

diname='di_record_inout';
fontsize=15;
if strcmp(diname,'di_record_inout')
    xxlim=[1e-10 10];
else
    xxlim=[1e-11 1];
end
yylim=0.1;
N=100;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:Nf
    
    di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname{isub},'/',fname{isub},'/',diname,'.txt']);

%% eliminate invalid value pour loglog diagram
    
%     if isub==1
%         Maxdi=max(abs(di_record_inout));
%     end
    di_norm=di_record_inout;%/Maxdi;
    % di_norm1=di_record_inout1/Maxdi;
    % di_norm2=di_record_inout2/Maxdi;
    
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
% 
%     %% Theo
%     plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/plel.txt']);
%     tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/tpel.txt']);
%     r2hill=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname{isub},'/r2hill_record.txt']);
% 
%     pla=plel(:,2);plamean=mean(pla);
%     ple=plel(:,3);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
%     plinc=plel(:,4);plincmean=mean(plinc);%plincstd=std(plinc);
%     
%     tpa=tpel(:,2);tpamean=mean(tpa);
%     tpe=tpel(:,3);tpemean=mean(tpe);
%     tpinc=tpel(:,4);tpincmean=mean(tpinc);
%     
%     N=1000000;%
%     %N=length(di_record_inout); %% sample size
%     %std0=abs(plincmean)^(1/3)+1;
%     if strcmp(diname,'de_record_inout')
%         [~,dinc]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
%     else
%         [dinc]=Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
%     end        
%     dincAbs=abs(dinc);
%     
%     %% compare MI
%     disp((sum(di_norm-mean(di_norm)).^2));
%     disp((sum(dinc-mean(dinc)).^2));
%     
    countx_p=histcounts(di_norm(di_norm>0),dix)'/length(di_norm(di_norm>0));
    countx_n=histcounts(-di_norm(di_norm<0),dix)'/length(-di_norm(di_norm<0));
%     countx_theo=histcounts(dincAbs,dix)'/length(dincAbs);
% 
    dixx=(dix(2:end)+dix(1:end-1))/2; 
    dixlen=dix(2:end)-dix(1:end-1);
    
    
%     x0=0.0000129761*0.5;
%     x0=1e-6;

%    x01=7.23133*10^-7;
%     f=@(x,x0)2*x0^2*x.*(x.^2 + x0^2).^(-2)/2;
%     f1=@(x,x0)abs(x).*(1 + x.^2).^(-2);    
%     plot_theo=f(dixx/180*pi,x01).*dixlen/180*pi;
%     plot_theo_p=plot_theo(dixx>=0);
%     plot_theo1=f1(dixx/180*pi/x01).*dixlen/x01/180*pi;
%     plot_theo_p1=plot_theo1(dixx>=0); 
    
    %countx=countx./dixlen;
    %% Direct distribution di DStb
%     aP=plamean;aT=tpamean;IP=plincmean;
%     eT=0;
%     mP=6.56e-9;
%     gm=3.5;
%     Rth=gm*aP*(mP/3)^(1/3);
%     
%     A=(3-aT/aP)/2;
%     B=(2-aT/aP)^(1/2);
%     dinc0=mP*aT/Rth*(2/(1-eT^2))^(1/2)/sqrt(A-B*cosd(IP));
% 
%     f=@(z)2./(2*(1+z.^2).^(3/2));
%     ploty_theotheo=f(dixx/180*pi/dinc0).*dixlen/dinc0/180*pi;


    plotx_p=dixx(dixx>=0);
    plotx=dixx;
    ploty_p=countx_p;
%     plotx_n=-dixx(dixx<0);
    ploty_n=countx_n;
    
%     ploty_theo=countx_theo;

    % ploty_p1=countx1(dix>=0);
    % ploty_n1=countx1(dix<0);
    % 
    % ploty_p2=countx2(dix>=0);
    % ploty_n2=countx2(dix<0);

    subplot(2,2,isub);
    
    semilogx(0,0,'w');hold all;
    %set(gca,'xscal','log');
%     bar(plotx_p,ploty_p,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
%     bar(plotx_n,ploty_n,'facecolor',[0.9 0.9 0.9],'edgecolor','none');

    pp=semilogx(plotx,ploty_p,'r.-','linewidth',1.3,'markersize',10);
    %bar(plotx_n,ploty_n,'facecolor',[0.9 0.9 0.9],'edgecolor','none');%axis square;
    pn=semilogx(plotx,ploty_n,'b.-','linewidth',1.3,'markersize',10);
%     semilogx(plotx,ploty_theo,'c-','linewidth',1.3,'markersize',10);
%     semilogx(plotx,ploty_theotheo,'y-','linewidth',1.3,'markersize',10);
% 
%     semilogx(plotx_p,plot_theo_p1,'g--','linewidth',1.3,'markersize',10);
%     semilogx(plotx_p,plot_theo_p,'g.-','linewidth',1.3,'markersize',10);

    % semilogx(plotx_p,ploty_p1,'r.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n1,'r.-','linewidth',1.3,'markersize',10);
    % 
    % semilogx(plotx_p,ploty_p2,'g.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n2,'g.-','linewidth',1.3,'markersize',10);


    %semilogy(plotx,ploty,'k.--','linewidth',1.3,'markersize',10);
    xlim(xxlim);
    ylim([0 yylim]);
    %ylim=get(gca,'ylim');

    % Ex=sum(plotx.*ploty);
    % %% Variance
    % Dx=sum((plotx-Ex).^2.*ploty);
    % Dx0=((plotx(end)-plotx(1))^2/12);
    % plot([Ex Ex],[ylim(1) ylim(2)],'k--');
    absMeandi=abs(mean(di_norm));
    Vardi=var(di_norm);
    plot([absMeandi absMeandi],[0 yylim],'k--');

    Meanabsdi=mean(abs(di_norm));
    Varabsdi=var(abs(di_norm));
    plot([Meanabsdi Meanabsdi],[0 yylim],'k--','linewidth',3.0);

    % TrimMeanabsdi=trimmean(abs(di_norm),1);
    % plot([TrimMeanabsdi TrimMeanabsdi],[ylim(1) ylim(2)],'k--','linewidth',1.3);

    % Rmsabsdi=rms(abs(di_norm));
    % plot([Rmsabsdi Rmsabsdi],[ylim(1) ylim(2)],'k--','linewidth',3.0);

%     Sumdi=abs(sum(di_norm));
%     plot([Sumdi Sumdi],[0 yylim],'k--','linewidth',3.0);

    hold off;
    if strcmp(diname,'di_record_inout')
        tag='I';
    else
        tag='e';
    end
    set(text(xxlim(1)*5,8/9*yylim,['$$N_{\Delta ',tag,'>0}= ',num2str(length(find(di_record_inout>=0))),'$$']),...,
        'Interpreter','latex','fontsize',fontsize,'color','r');
    set(text(xxlim(1)*5,7/9*yylim,['$$N_{\Delta ',tag,'<0}= ',num2str(length(find(di_record_inout<0))),'$$']),...,
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
