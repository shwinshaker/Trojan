%%bar di
clear;

ffname1='RealPlutinosNpl';
ffname2='RanPlutinosNpl';
fname1={'1999CE119';'2001FU172';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};%;'2001FU172';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};%;'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
suf='40pl';
dd='di';

fontsize=15;
yylim=0.2;
N=50;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:length(fname1)
%     if isub==3
%       suf='40pl'; 
%     else
%       suf='20pl';
%     end
if strcmp(dd,'di')   
    di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',[fname1{isub},'_1Gyr_',suf],'/di_record_inout.txt'));
    di_record_inout2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',[fname1{isub},'_1Gyr_',suf,'_ran'],'/di_fit_perturb.txt'));
    %di_record_inout2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',[fname1{isub},'_40plTimes'],'/di_fit_perturb.txt'));
else 
    di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',[fname1{isub},'_1Gyr_',suf],'/de_record_inout.txt'));
    di_record_inout2=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname2,'/',[fname1{isub},'_1Gyr_',suf,'_ran'],'/de_fit_perturb.txt'));
end


%% eliminate invalid value pour loglog diagram
    
%     if isub==1
%         Maxdi=max(abs(di_record_inout));
%     end
    di_norm=di_record_inout;%/Maxdi;
    di_norm2=di_record_inout2;
    % di_norm2=di_record_inout2/Maxdi;
 
    di_norm_0=di_norm(di_norm~=0);
    di_norm_02=di_norm2(di_norm2~=0);
    
    H=kstest2(di_norm,di_norm2);
    disp(titlename{isub});
    disp(H);
    
    %if isub==1
        Minlog=-log(min(min(abs(di_norm_0)),min(abs(di_norm_02))));

        dix=0:N;
        dix=dix';
        dlog=Minlog/N;
        dix=exp(-dix*dlog);
        %dix=[-dix;dix];
        dix=sort(dix);
    %end

    %countx=histcounts(di_norm,dix)'/length(di_record_inout);

    dixx=(dix(2:end)+dix(1:end-1))/2;

%     plotx_p=dixx(dixx>=0);
%     ploty_p=countx(dixx>=0);
%     plotx_n=-dixx(dixx<0);
%     ploty_n=countx(dixx<0);

    countx_real=histcounts(abs(di_norm_0),dix)'/length(di_record_inout);
    countx_ran=histcounts(abs(di_norm_02),dix)'/length(di_record_inout2);

    % ploty_p1=countx1(dix>=0);
    % ploty_n1=countx1(dix<0);
    % 
    % ploty_p2=countx2(dix>=0);
    % ploty_n2=countx2(dix<0);

    subplot(2,2,isub);
    
    semilogx(0,0,'w');hold all;
    %set(gca,'xscal','log');
    
    h1=semilogx(dixx,countx_real,'k.-','linewidth',1.3,'markersize',10);
    h2=semilogx(dixx,countx_ran,'r.-','linewidth',1.3,'markersize',10);
    legend([h1 h2],{'$Integ$','$Ran$'},'box','off','Interpreter','latex','fontsize',fontsize/3*2);

%     semilogx(plotx_p,ploty_p,'r.-','linewidth',1.3,'markersize',10);
%     semilogx(plotx_n,ploty_n,'b.-','linewidth',1.3,'markersize',10);

    % semilogx(plotx_p,ploty_p1,'r.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n1,'r.-','linewidth',1.3,'markersize',10);
    % 
    % semilogx(plotx_p,ploty_p2,'g.-','linewidth',1.3,'markersize',10);
    % semilogx(plotx_n,ploty_n2,'g.-','linewidth',1.3,'markersize',10);


    %semilogy(plotx,ploty,'k.--','linewidth',1.3,'markersize',10);
    xlim([1e-10 10]);
    ylim([0 yylim]);
    %ylim=get(gca,'ylim');

    % Ex=sum(plotx.*ploty);
    % %% Variance
    % Dx=sum((plotx-Ex).^2.*ploty);
    % Dx0=((plotx(end)-plotx(1))^2/12);
    % plot([Ex Ex],[ylim(1) ylim(2)],'k--');
    
%     absMeandi=abs(mean(di_norm));
%     Vardi=var(di_norm);
%     plot([absMeandi absMeandi],[0 yylim],'k--');

    Meanabsdi=mean(abs(di_norm));
    Varabsdi=var(abs(di_norm));
    plot([Meanabsdi Meanabsdi],[0 yylim],'k--','linewidth',1.0);
    
%     absMeandi2=abs(mean(di_norm2));
%     Vardi2=var(di_norm2);
%     plot([absMeandi2 absMeandi2],[0 yylim],'r--');

    Meanabsdi2=mean(abs(di_norm2));
    Varabsdi2=var(abs(di_norm2));
    plot([Meanabsdi2 Meanabsdi2],[0 yylim],'r--','linewidth',1.0);

    % TrimMeanabsdi=trimmean(abs(di_norm),1);
    % plot([TrimMeanabsdi TrimMeanabsdi],[ylim(1) ylim(2)],'k--','linewidth',1.3);

    % Rmsabsdi=rms(abs(di_norm));
    % plot([Rmsabsdi Rmsabsdi],[ylim(1) ylim(2)],'k--','linewidth',3.0);

    Sumdi=abs(sum(di_norm));
    plot([Sumdi Sumdi],[0 yylim],'k--','linewidth',3.0);
    
    Sumdi2=abs(sum(di_norm2));
    plot([Sumdi2 Sumdi2],[0 yylim],'r--','linewidth',3.0);
    
    
    hold off;

    set(text(5e-10,8/9*yylim,['$$N= ',num2str(length(di_record_inout)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');
    set(text(5e-10,7/9*yylim,['$$H_{KS}= ',num2str(H),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');

    %set(text(5e-10,7/9*yylim,['$$N_{di<0}= ',num2str(length(find(di_record_inout<0))),'$$']),'Interpreter','latex','fontsize',fontsize,'color','b');
    title(titlename{isub},'fontsize',fontsize);

    set(gca,'xTick',power(10,-10:1:2));
    
    if strcmp(dd,'di')   
        xlabel('$\Delta i~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    else 
        xlabel('$\Delta e$','fontsize',fontsize,'Interpreter','latex');
    end
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end
