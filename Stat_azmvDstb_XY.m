clear;
clear;

ffname1='RealPlutinosNpl';
fname1={'1999CE119';'2001FU172';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

suf='40pl';

fontsize=15;
N=300;
yylim=0.1;

figure(1);
set(gcf,'Position',[400,100,700,500],'color','w');

for isub=1:4
    subplot(2,2,isub);
    plot(0,0,'w');hold all;
        
    for iplot=1:1
        ffname=eval(['ffname',num2str(iplot)]);
        if iplot==1
            fname=[fname1{isub},'_1Gyr_',suf];
            titlename=titlename1{isub};
            color='k';
        else
            %fname=[fname1{isub},'_40plTimes']; 
            fname=[fname1{isub},'_1Gyr_',suf,'_ran'];
            color='r';
        end
        
    r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
    r2hill_mean=mean(r2hill_record);
    hill=sqrt(r2hill_mean);
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));
    AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_pl.txt'));
    AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_tp.txt'));

    %di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/di_record_inout.txt'));

    Vp=PV_record_pl(:,4:6);
    Vt=PV_record_tp(:,4:6);
    PV_rlt=PV_record_pl-PV_record_tp;
    P_rlt=PV_rlt(:,1:3);
    V_rlt=PV_rlt(:,4:6);
    Vrnorm=(V_rlt(:,1).^2+V_rlt(:,2).^2+V_rlt(:,3).^2).^(1/2);
    Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
    sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
    coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));
    P_rlt_con=zeros(length(P_rlt),3);
    V_rlt_con=zeros(length(V_rlt),3);
    P_rlt_con_norm=zeros(length(P_rlt),1);
    Vp_con=zeros(length(P_rlt),3);
    Vt_con=zeros(length(P_rlt),3);
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
         Vp_con(i,:)=(P3*P2*P1*Vp(i,:)')';
         Vt_con(i,:)=(P3*P2*P1*Vt(i,:)')';
    end
    
    xv=(V_rlt_con(:,1).^2+V_rlt_con(:,2).^2).^(1/2);
    yv=V_rlt_con(:,3);
    thetavz_rlt=atan(yv./xv);
%     thetavz_rlt=atan(V_rlt_con(:,2)./V_rlt_con(:,1));
    
%     xv=(Vp_con(:,1).^2+Vp_con(:,2).^2).^(1/2);
%     yv=Vp_con(:,3);
%     thetavp=atan(yv./xv);
%     thetavp=atan(Vp_con(:,2)./Vp_con(:,1));

    thetax=0:N;
    thetax=thetax';
    dtheta=(pi/2--pi/2)/N;
    thetax=-pi/2+thetax*dtheta;
    
    countvz_rlt=histcounts(thetavz_rlt,thetax)';
%     countvp=histcounts(thetavp,thetax)';

    plotx=(thetax(2:end)-dtheta/2)/pi*180;
    plotyvz_rlt=countvz_rlt/length(thetavz_rlt);
%     plotyvp=countvp/length(thetavp);
    
    disp(titlename);

%     bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
%     h2=plot(plotx,plotyvp,'k.-','linewidth',1.3,'markersize',10);
    h1=plot(plotx,plotyvz_rlt,'k.-','linewidth',1.3,'markersize',10);

    xlim([-90 90]);
    ylim([0 yylim]);

    %% Expectation
    %Ex=sum(plotx.*ploty);
    %% Variance
    %Dx=sum((plotx-Ex).^2.*ploty);
    %Dx0=((plotx(end)-plotx(1))^2/12);
    %plot([Ex Ex],[ylim(1) ylim(2)],'k--');

%     %%Mean
%     MeanAzm=mean(thetav)/pi*180;
%     VarAzm=var(thetav);
%     plot([MeanAzm MeanAzm],[0 yylim],[color,'--']);
    
    end
    hold off;
%     legend([h1 h2],{'$v_r$','$v_p$'},'box','off','Interpreter','latex','fontsize',fontsize/3*2);

        
    %set(text(20,8/9*yylim,['$$N_{CE}= ',num2str(length(theta)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     set(text(-20,8/9*yylim,['$$\overline{\theta_v}= ',num2str(MeanAzm,4),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
    title(titlename,'fontsize',fontsize);

    set(gca,'xTick',-90:30:90);
    xlabel('$\theta_v~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end


