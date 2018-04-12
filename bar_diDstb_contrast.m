%%bar di
clear;
% ffname='RealPlutinosNpl';
% %ffname2='RealPlutinosNpl';
% ffname1='RanPlutinos';
% ffname2='RealPlutinos';
% fname1='1999CE119_10k';
% fname='1999CE119_1Gyr_40pl';
% %fname2='2001KN77_1Gyr_40pl';
% fname2='1999CE119_1Gyr';
% %fname2='2001FU172_1Gyr_40pl';

ffname='RealPlutinosNpl';
fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

fontsize=15;
yylim=0.1;
N=50;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    
%PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
%PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));

    di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/di_record_inout.txt'));
% di_record_inout2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',fname2,'/di_record_inout.txt'));
% %di_record_inout0=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname0,'/',fname0,'/di_record_inout.txt'));
% 
% %di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/di_fit_perturb.txt'));
% di_record_inout1=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname1,'/di_fit_perturb.txt'));
% 
%PV_rlt=PV_record_tp-PV_record_pl;
%r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
%r2hill_mean=mean(r2hill_record);
%hill=sqrt(r2hill_mean);



%% eliminate invalid value pour loglog diagram

    Maxdi=max(abs(di_record_inout));
    di_norm=di_record_inout/Maxdi;
    % di_norm1=di_record_inout1/Maxdi;
    % di_norm2=di_record_inout2/Maxdi;

    
    di_norm_0=di_norm(di_norm~=0);
    
    H=kstest2(di_norm_0(di_norm_0>0),-di_norm_0(di_norm_0<0));
    disp(titlename{isub});
    disp(H);
    
    Minlog=-log(min(abs(di_norm_0)));

    dix=0:N;
    dix=dix';
    dlog=Minlog/N;
    dix=exp(-dix*dlog);
    dix=[-dix;dix];
    dix=sort(dix);

    countx=histcounts(di_norm,dix)'/length(di_record_inout);

    % countx1=zeros(length(dix)-1,1);
    % for i=1:length(dix)-1
    %     countx1(i)=length(find(dix(i) < di_norm1 & di_norm1 <= dix(i+1)));
    % end
    % countx1=countx1/length(di_record_inout1);
    % 
    % countx2=zeros(length(dix)-1,1);
    % for i=1:length(dix)-1
    %     countx2(i)=length(find(dix(i) < di_norm2 & di_norm2 <= dix(i+1)));
    % end
    % countx2=countx2/length(di_record_inout2);

    dix=(dix(2:end)+dix(1:end-1))/2;

    plotx_p=dix(dix>=0);
    ploty_p=countx(dix>=0);
    plotx_n=-dix(dix<0);
    ploty_n=countx(dix<0);

    % ploty_p1=countx1(dix>=0);
    % ploty_n1=countx1(dix<0);
    % 
    % ploty_p2=countx2(dix>=0);
    % ploty_n2=countx2(dix<0);

    subplot(2,2,isub);
    
    semilogx(0,0,'w');hold all;
    %bar(plotx_p,ploty_p,'facecolor',[0.9 0.9 0.9],'edgecolor','none');%axis square;
    %set(gca,'xscal','log');
    semilogx(plotx_p,ploty_p,'r.-','linewidth',1.3,'markersize',10);
    %bar(plotx_n,ploty_n,'facecolor',[0.9 0.9 0.9],'edgecolor','none');%axis square;
    semilogx(plotx_n,ploty_n,'b.-','linewidth',1.3,'markersize',10);

    % semilogx(plotx_p,ploty_p1,'r.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n1,'r.-','linewidth',1.3,'markersize',10);
    % 
    % semilogx(plotx_p,ploty_p2,'g.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n2,'g.-','linewidth',1.3,'markersize',10);


    %semilogy(plotx,ploty,'k.--','linewidth',1.3,'markersize',10);
    xlim([1e-8 10]);
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
    plot([Meanabsdi Meanabsdi],[0 yylim],'k--','linewidth',2.0);

    % TrimMeanabsdi=trimmean(abs(di_norm),1);
    % plot([TrimMeanabsdi TrimMeanabsdi],[ylim(1) ylim(2)],'k--','linewidth',1.3);

    % Rmsabsdi=rms(abs(di_norm));
    % plot([Rmsabsdi Rmsabsdi],[ylim(1) ylim(2)],'k--','linewidth',3.0);

%     Sumdi=abs(sum(di_norm));
%     plot([Sumdi Sumdi],[0 yylim],'k--','linewidth',3.0);

    hold off;

    set(text(2e-8,8/9*yylim,['$$N_{di>0}= ',num2str(length(find(di_record_inout>=0))),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');
    set(text(2e-8,7/9*yylim,['$$N_{di<0}= ',num2str(length(find(di_record_inout<0))),'$$']),'Interpreter','latex','fontsize',fontsize,'color','b');
    title(titlename{isub},'fontsize',fontsize);

    set(gca,'xTick',power(10,-8:1:2));

    xlabel('$\Delta i~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end
