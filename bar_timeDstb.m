%% time distribution bar
%% ran&record contrast
clear;
ffname='RealPlutinos';
fname={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

fontsize=15;
N=50;
yylim1=0.1;
yylim2=1.0;

set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    subplot(2,2,isub);

    CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/CE_record.txt'));
    %CE_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ran_record.txt'));

    record_time=CE_record(:,1);
    %maxtime=max(el(:,1));
    maxtime=max(record_time);
    %ran_time=ran_record(:,1)/ran_record(end,1)*maxtime; %% time scale
    time=record_time/maxtime;
    
    [a,b]=unifit(time);
    cdf=unifcdf(time,a,b);
%     cdf=unifrnd(0,1,length(time),1);
%     cdf=sort(cdf);
    [H,p,ksstat,cv]=kstest(time,[time,cdf]);
    disp(titlename{isub});
    disp([H p ksstat cv]);
    disp(kstest2(time,cdf));
    
    disp(min(CE_record(:,2)));
    
    [f,x_value]=ecdf(time);
    
    timebar=0:N;
    timebar=timebar';
    timebar=timebar/N;
    
    plotx=timebar(2:end)-1/2/N;
    
    count=histcounts(time,timebar)';
    
    ploty=count/length(record_time);
    
    [hax,hline1,hline2]=plotyy(plotx,ploty,[x_value,x_value],[f,unifcdf(x_value,a,b)]);hold all;
    hline1.Color='w';
    hline2(1).Color='k';
    hline2(2).Color='r';

    bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');%axis square;
    plot(plotx,ploty,'k.--','linewidth',1.3,'markersize',10);
%     plot(x_value,f,'r.--','linewidth',1.0,'markersize',8);
%     plot(x_value,unifcdf(x_value,a,b),'k.--','linewidth',1.0,'markersize',8);


    plot([mean(time) mean(time)],[0 yylim2],'k--');
    
    ylim(hax(1),[0,yylim1]);
    ylim(hax(2),[0,yylim2]);
    set(hax(1),'ycolor','k');
    set(hax(1),'ytick',0:yylim1/5:yylim1);
    set(hax(2),'ytick',0:0.2:1);

%     ylim([0 yylim]);
%     xlim([0 1]);
    
    hold off;

    set(text(0.1,8/9*yylim1,['$$N_{CE}= ',num2str(length(record_time)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename{isub},'fontsize',fontsize);
    xlabel('$\mathrm{t/t_{tot}}$','fontsize',fontsize,'Interpreter','latex')
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex')
end
%axis square;

%axis([0 N 0 max(count)]);
% set(gca,'xticklabel',0:estep:N*estep);
%set(text(0.5,8/9*(max(count(x))),strcat('N_{CE}= ',num2str(length(CE_record)))),'color','red');
%title(num2str(length(CErecord),'N=%04d'));
%title('CE-bar');
% filename=num2str(num,'20141%03 d.fits'); 

% subplot(2,2,4);
% plot(x/N,countran(x),'k-');
% %axis square;
% set(gca,'ylim',ylim_1);
% % set(text(0.5,8/9*(max(count(x))),strcat('N_{CE}= ',num2str(length(ran_CE_di)))),'color','red');
% % title('ran-CE-bar');