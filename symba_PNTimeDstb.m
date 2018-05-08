%% time distribution bar
%% ran&record contrast
clear;
ffname='symbaRealPlutinosNpl_fast';
fname={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

diname='di_fit_perturb';

fontsize=15;
N=50;
yylim1=0.05;
yylim2=1.0;

set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    subplot(2,2,isub);

    CE_record=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/CE_record.txt']);
    [C,ia,ic]=unique(CE_record(:,3));
    CE_record=CE_record(ia,:);
    di=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',diname,'.txt']);

    record_time=CE_record(:,1);
    %maxtime=max(el(:,1));
    maxtime=max(record_time);
    %ran_time=ran_record(:,1)/ran_record(end,1)*maxtime; %% time scale
    time=record_time/maxtime;
    
%     [a,b]=unifit(time);
%     cdf=unifcdf(time,a,b);
% %     cdf=unifrnd(0,1,length(time),1);
% %     cdf=sort(cdf);
%     [H,p,ksstat,cv]=kstest(time,[time,cdf]);
%     disp(titlename{isub});
%     disp([H p ksstat cv]);
%     disp(kstest2(time,cdf));
%     
%     disp(min(CE_record(:,2)));
    

    
    timebar=0:N;
    timebar=timebar';
    timebar=timebar/N;
    
    plotx=timebar(2:end)-1/2/N;
    
    time_p=time(di>0);
    time_n=time(di<0);
    [f_p,x_value_p]=ecdf(time_p);
    [f_n,x_value_n]=ecdf(time_n);
        
    count_p=histcounts(time_p,timebar)';
    count_n=histcounts(time_n,timebar)';
    
    ploty_p=count_p/length(record_time);
    ploty_n=count_n/length(record_time);
    
    
%     [hax,hline1,hline2]=plotyy([plotx,plotx],[ploty_p,ploty_n],[x_value_p;NaN;x_value_n],[f_p;NaN;f_n]);hold all;
    [hax1,hline11,hline21]=plotyy(plotx,ploty_p,x_value_p,f_p);hold all;
    [hax2,hline12,hline22]=plotyy(plotx,ploty_n,x_value_n,f_n);

    hline11.Color='w';
    hline12.Color='w';
    hline21.Color='r';
    hline22.Color='b';

    bar(plotx,ploty_p,'facecolor',[0.95 0.95 0.95],'edgecolor','none');%axis square;
    bar(plotx,ploty_n,'facecolor',[0.95 0.95 0.95],'edgecolor','none');%axis square;
    
    plot(plotx,ploty_p,'r.-','linewidth',1.2,'markersize',10);
    plot(plotx,ploty_n,'b.-','linewidth',1.2,'markersize',10);

    plot([mean(time_p) mean(time_p)],[0 yylim2],'r--');
    plot([mean(time_n) mean(time_n)],[0 yylim2],'b--');

    ylim(hax1(1),[0,yylim1]);
    ylim(hax1(2),[0,yylim2]);
    set(hax1(1),'ycolor','k');
    set(hax1(2),'ycolor','k');
    set(hax2(1),'ycolor','k');
    set(hax2(2),'ycolor','k');
    set(hax1(1),'ytick',[]);
    set(hax1(2),'ytick',[]);
    set(hax2(1),'ytick',0:yylim1/5:yylim1);
    set(hax2(2),'ytick',0:0.2:1);

%     ylim([0 yylim]);
%     xlim([0 1]);
    
    hold off;

    set(text(0.1,8/9*yylim1,['$$N_{CE}= ',num2str(length(di)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename{isub},'fontsize',fontsize);
    xlabel('$t/t_{tot}$','fontsize',fontsize,'Interpreter','latex')
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