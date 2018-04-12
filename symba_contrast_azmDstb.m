clear;
clear;
ffname1='RealPlutinosNpl';
ffname2='symbaRealPlutinosNpl_fast';
fname1={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
fname2={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

fontsize=15;
Nbin=50;
yylim=0.04;
% xxlim=180;
xxlim1=-90;xxlim2=90;

set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:length(fname1)
    subplot(2,2,isub);
    plot(0,0,'w');hold all;
    
    for iplot=1:2
        switch iplot
            case 1
                color='k';
            case 2
                color='r';
        end
        ffname=eval(['ffname',num2str(iplot)]);
        disp(ffname);
        fname=eval(['fname',num2str(iplot)]);
%         r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
%         r2hill_mean=mean(r2hill_record);
%         hill=sqrt(r2hill_mean);
        PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_pl.txt'));
        PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_tp.txt'));
        AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_pl.txt'));
        AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_tp.txt'));
        
        %di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/di_record_inout.txt'));
        
        PV_rlt=PV_record_pl-PV_record_tp;
        P_rlt=PV_rlt(:,1:3);
        V_rlt=PV_rlt(:,4:6);
        Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
        sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
        coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));
        P_rlt_con=zeros(length(P_rlt),3);
        V_rlt_con=zeros(length(V_rlt),3);
        P_rlt_con_norm=zeros(length(P_rlt),1);
        for i=1:length(P_rlt)
            P1=[cosd(AE_record_tp(i,4)) sind(AE_record_tp(i,4)) 0; ...,
                -sind(AE_record_tp(i,4)) cosd(AE_record_tp(i,4)) 0; ...,
                0 0 1];
            
            P2=[1 0 0;...,
                0 cosd(AE_record_tp(i,3)) sind(AE_record_tp(i,3)); ...,
                0 -sind(AE_record_tp(i,3)) cosd(AE_record_tp(i,3))];
            
            P3=[coswf(i) sinwf(i) 0;...,
                -sinwf(i) coswf(i) 0;...,
                0 0 1];
            P_rlt_con(i,:)=(P3*P2*P1*P_rlt(i,:)')';
            V_rlt_con(i,:)=(P3*P2*P1*V_rlt(i,:)')';
            P_rlt_con_norm(i)=norm(P_rlt_con(i,:));
        end
        clear i;
        
        x=(P_rlt_con(:,1).^2+P_rlt_con(:,2).^2).^(1/2);
        y=P_rlt_con(:,3);
        
%         x=P_rlt_con(:,1);
%         y=P_rlt_con(:,2);
        % check i
        if isreal(i)
            error('Fatal error! i is assigned!');
        end
        ang=circ_rad2ang(angle(x+y*i));
        
        eval(['theta',num2str(iplot),'=ang;']);
        theta=eval(['theta',num2str(iplot)]);
%         for ii=1:length(theta)
%             if x(ii)>0 & y(ii)>0
%                 theta(ii)=theta(ii);
%             elseif x(ii)<0 & y(ii)>0
%                 theta(ii)=pi+theta(ii);
%             elseif x(ii)<0 & y(ii)<0
%                 theta(ii)=theta(ii)-pi;
%             elseif x(ii)>0 & y(ii)<0
%                 theta(ii)=theta(ii);
%             end
%         end
        
        
        disp(titlename{isub});
        
        [binx,plotx,dx]=Fun_bin(xxlim1,xxlim2,Nbin);
        
%         thetax=0:N;
%         thetax=thetax';
%         dtheta=(pi--pi)/N;
%         thetax=-pi+thetax*dtheta;
%         
        countx=histcounts(theta,binx)'/dx;
%         plotx=(thetax(2:end)-dtheta/2)/pi*180;
        ploty=countx/length(theta);
        
        %bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        h=plot(plotx,ploty,[color,'.-'],'linewidth',1.3,'markersize',10);
        eval(['h',num2str(iplot),'=h;']);
        
        %% stamp the peak
        %     if isub==1
        %         [peak,peaki]=findpeaks(ploty,'MinPeakHeight',0.02);
        %         %plot(plotx(peaki),ploty(peaki),'r.','linewidth',1.3,'markersize',10);
        %         for ii=1:length(peaki)
        %             ipeak=peaki(ii);
        %             if ploty(ipeak)-ploty(ipeak+1)<ploty(ipeak)-ploty(ipeak-1)
        %                 inear=ipeak+1;
        %             else
        %                 inear=ipeak-1;
        %             end
        %             pplotx=(plotx(inear)+plotx(ipeak))/2;
        %             plot([pplotx pplotx],[0 yylim],'r-','linewidth',1.3);
        %         end
        %     end
        
        xlim([xxlim1 xxlim2]);
        ylim([0 yylim]);
        
        %% Expectation
        %Ex=sum(plotx.*ploty);
        %% Variance
        %Dx=sum((plotx-Ex).^2.*ploty);
        %Dx0=((plotx(end)-plotx(1))^2/12);
        %plot([Ex Ex],[ylim(1) ylim(2)],'k--');
        
        %     %%Mean
        %     MeanAzm=mean(theta)/pi*180;
        %     VarAzm=var(theta);
        %     plot([MeanAzm MeanAzm],[0 yylim],[color,'--']);
        
        if iplot==1
            plot([90 90],[0 yylim],[color,'--']);
            plot([-90 -90],[0 yylim],[color,'--']);
            plot([0 0],[0 yylim],[color,'--']);
        end
        
    end
    
    hold off;
    
    %     H=kstest2(theta1,theta2);
    %     disp(H);
    legend([h1 h2],{'Rmvs','Symba'},'box','off','Interpreter','latex','fontsize',fontsize/5*4);
    
    %set(text(-75,8/9*yylim,['$$H_{KS}= ',num2str(H),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');
    
    %set(text(20,8/9*yylim,['$$N_{CE}= ',num2str(length(theta)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    %set(text(-20,8/9*yylim,['$$\overline{\theta}= ',num2str(MeanAzm,4),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename{isub},'fontsize',fontsize);
    
    set(gca,'xTick',xxlim1:30:xxlim2);
    xlabel('$\varphi_0~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end




