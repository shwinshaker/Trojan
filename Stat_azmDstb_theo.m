
if ~exist('theta_dat1','var')
    
    clear;
    %% data
    ffname='RealPlutinosNpl';
    ffname_theo='RealPlutinos';
    fnameL={'1999CE119';'2001FU172';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
    titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
    suf='40pl';
    
    inda=2;
    
    for ipl=1:4
        fname=[fnameL{ipl},'_1Gyr_',suf];
        titlename=titlename1{ipl};
        
        disp(titlename);
        PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
        PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));
        AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_pl.txt'));
        AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_tp.txt'));
        
        ind=find(AE_record_pl(:,1)>39 & AE_record_pl(:,1)<40);
        PV_record_pl=PV_record_pl(ind,:);
        PV_record_tp=PV_record_tp(ind,:);
        AE_record_pl=AE_record_pl(ind,:);
        AE_record_tp=AE_record_tp(ind,:);
        
        Vp=PV_record_pl(:,4:6);
        Vt=PV_record_tp(:,4:6);
        PV_rlt=PV_record_pl-PV_record_tp;
        P_rlt=PV_rlt(:,1:3);
        V_rlt=PV_rlt(:,4:6);
        Vrnorm=(V_rlt(:,1).^2+V_rlt(:,2).^2+V_rlt(:,3).^2).^(1/2);
        Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
        sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
        coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));
        cosf=coswf.*cosd(AE_record_tp(:,5))+sinwf.*sind(AE_record_tp(:,5));
        sinf=sinwf.*cosd(AE_record_tp(:,5))-coswf.*sind(AE_record_tp(:,5));
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
        
        %% Theta_P
        xv=(Vp_con(:,1).^2+Vp_con(:,2).^2).^(1/2);
        yv=Vp_con(:,3);
        thetap_dat=atand(yv./xv);
        eval(['thetap_dat',num2str(ipl),'=thetap_dat;']);
        
        %% Phi_P
        phip_dat=atand(Vp_con(:,2)./Vp_con(:,1));
        x=Vp_con(:,1);
        y=Vp_con(:,2);
        for ii=1:length(phip_dat)
            if x(ii)>0 && y(ii)>0
                phip_dat(ii)=phip_dat(ii);
            elseif x(ii)<0 && y(ii)>0
                phip_dat(ii)=180+phip_dat(ii);
            elseif x(ii)<0 && y(ii)<0
                phip_dat(ii)=phip_dat(ii)-180;
            elseif x(ii)>0 && y(ii)<0
                phip_dat(ii)=phip_dat(ii);
            end
        end
        eval(['phip_dat',num2str(ipl),'=phip_dat;']);
        
        
        %% phi_vPvT
        phivpvt_dat=acosd(sum(Vp_con.*Vt_con,2)./(sum(Vp_con.*Vp_con,2).^(1/2).*sum(Vt_con.*Vt_con,2).^(1/2)));
        eval(['phivpvt_dat',num2str(ipl),'=phivpvt_dat;']);
        
        %% Theta_v
        xv=(V_rlt_con(:,1).^2+V_rlt_con(:,2).^2).^(1/2);
        yv=V_rlt_con(:,3);
        thetav_dat=atand(yv./xv);
        xv=V_rlt_con(:,1);
        yv=V_rlt_con(:,2);
        %     figure;
        %     plot(xv,yv,'k.');
        phiv_dat=atand(yv./xv);
        phiv_dat(xv<0 & yv<0)=phiv_dat(xv<0 & yv<0)+180;
        phiv_dat(xv<0 & yv>0)=phiv_dat(xv<0 & yv>0)+180;
        %     phiv_dat(xv>0 & yv<0)=phiv_dat(xv>0 & yv<0)+360;
        
        eval(['thetav_dat',num2str(ipl),'=thetav_dat;']);
        eval(['phiv_dat',num2str(ipl),'=phiv_dat;']);
        
        %% Theta
        xv=(P_rlt_con(:,1).^2+P_rlt_con(:,2).^2).^(1/2);
        yv=P_rlt_con(:,3);
        theta_dat=atand(yv./xv);
        eval(['theta_dat',num2str(ipl),'=theta_dat;']);
        xv=P_rlt_con(:,1);
        yv=P_rlt_con(:,2);
        phi_dat=atand(yv./xv);
        phi_dat(xv<0 & yv<0)=phi_dat(xv<0 & yv<0)+180;
        phi_dat(xv<0 & yv>0)=phi_dat(xv<0 & yv>0)+180;
        phi_dat(xv>0 & yv<0)=phi_dat(xv>0 & yv<0)+360;
        eval(['phi_dat',num2str(ipl),'=phi_dat;']);
        
        
        %% lambda
        coswf(coswf>1)=1;
        coswf(coswf<-1)=-1;
        lambda_dat=acosd(coswf);
        indqd4=coswf>0 & sinwf<0;
        lambda_dat(indqd4)=360-lambda_dat(indqd4);
        indqd3=coswf<0 & sinwf<0;
        lambda_dat(indqd3)=360-lambda_dat(indqd3);
        eval(['lambda_dat',num2str(ipl),'=lambda_dat;']);
        
        %% f
        cosf(cosf>1)=1;
        cosf(cosf<-1)=-1;
        f_dat=acosd(cosf);
        ind=cosf>0 & sinf<0;
        f_dat(ind)=360-f_dat(ind);
        ind=cosf<0 & sinf<0;
        f_dat(ind)=360-f_dat(ind);
        eval(['f_dat',num2str(ipl),'=f_dat;']);
        
    end

end

if ~exist('theta1','var')
    
    for ipl=1:4
        %% Theo theta Distribution
        %% get plel
        plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname_theo,'/',fnameL{ipl},'_1Gyr/plel.txt']);
        tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname_theo,'/',fnameL{ipl},'_1Gyr/tpel.txt']);
    
        pla=plel(:,inda);plamean=mean(pla);
        ple=plel(:,inda+1);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
        plinc=plel(:,inda+2);plincmean=mean(plinc);%plincstd=std(plinc);
        
        tpa=tpel(:,inda);tpamean=mean(tpa);
        tpe=tpel(:,inda+1);tpemean=mean(tpe);
        tpinc=tpel(:,inda+2);tpincmean=mean(tpinc);
        
        N=1000000;
        %% ?????????
        %     std0=abs(plincmean-tpincmean)^(1/3)+1;
        % %     disp(std0);
        %     [~,~,thetap,phip,thetav,phiv,phivpvt,lambda,theta,phi,f] = ...,
        %         Fun_diDstb_theo(N,1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
        std0=abs(plincmean)^(1/3)+1;
        [theta,phi,thetap,phip,thetav,phiv,phivpvt,lambda,f]=...,
            Fun_azmDstb_theo(N,plamean,plemean,plincmean,tpamean,tpincmean,std0);
        
        eval(['thetap',num2str(ipl),'=thetap;']);
        eval(['phip',num2str(ipl),'=phip;']);
        eval(['phivpvt',num2str(ipl),'=phivpvt;']);
        eval(['thetav',num2str(ipl),'=thetav;']);
        eval(['phiv',num2str(ipl),'=phiv;']);
        eval(['theta',num2str(ipl),'=theta;']);
        eval(['phi',num2str(ipl),'=phi;']);
        eval(['lambda',num2str(ipl),'=lambda;']);
        eval(['f',num2str(ipl),'=f;']);
        
    end

end

%% bin and plot
figure;
set(gcf,'Position',[400,100,700,500],'color','w');

Nbin=50; 
Nbin_theo=400;
%% default
xxlim1=-90;xxlim2=90;
yylim=0.2;
dxxlim=30;
fontsize=15;

varname='theta';%'thetap','phip','phivpvt','thetav','lambda'

%% lambda ???????????0?180????
switch varname
    case 'thetap'
        xxlim1=-60;xxlim2=60;
        label='\theta_{v_P}';  
        yylim=0.15;
    case 'phip'
        xxlim1=30;xxlim2=150;
        yylim=0.1;
        label='\varphi_{v_P}';
    case 'phivpvt'
        xxlim1=0;xxlim2=90;
        yylim=0.2;
        label='\varphi_v';
    case 'thetav'
        label='\theta_{v_r}';
        yylim=0.04;
    case 'phiv'
        label='\varphi_{v_r}';
        yylim=0.04;
        xxlim1=-90;xxlim2=270;
    case 'theta'
        yylim=0.04;
        label='\theta_0';
    case 'phi'
        label='\varphi_0';
        yylim=0.02;
        xxlim1=0;xxlim2=360;
    case 'f'
        label='f';
        yylim=0.02;
        xxlim1=0;xxlim2=360;
    case 'lambda'
        xxlim1=0;xxlim2=360;
        dxxlim=60;
        yylim=0.01;
        label='\widetilde{\lambda_T}';
%     otherwise
end

[countx,plotx,dx]=Fun_bin(xxlim1,xxlim2,Nbin);
[countx_theo,plotx_theo,dx_theo]=Fun_bin(xxlim1,xxlim2,Nbin_theo);

for isub=1:4
    subplot(2,2,isub)
    var_teo=eval([varname,num2str(isub)]);
    var_dat=eval([varname,'_dat',num2str(isub)]);

    %% data
    counts=histcounts(var_dat,countx)'/dx;
    ploty=counts/length(var_dat);
    plot(0,0,'w+');hold all;
    bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    h1=plot(plotx,ploty,'k.--','linewidth',1.3,'markersize',10);hold all;
    
    %% theo
    counts=histcounts(var_teo,countx_theo)'/dx_theo;
    ploty=counts/length(var_teo);
    h2=plot(plotx_theo,ploty,'r-','linewidth',1.3);
    
%     h2=plot([xxlim1 xxlim2],[1/360 1/360],'r-','linewidth',1.3);

    plot([90 90],[0 yylim],'k--');
    plot([270 270],[0 yylim],'k--');
    plot([0 0],[0 yylim],'k--');
    plot([180 180],[0 yylim],'k--');
    
    xlim([xxlim1 xxlim2]);
    ylim([0 yylim]);
    set(gca,'xTick',xxlim1:dxxlim:xxlim2);
    
%     plot([max(plinc) max(plinc)],[0 1],'k--');
%     plot([-max(plinc) -max(plinc)],[0 1],'k--');
    xlabel(['$',label,'~\mathrm{(DEG)}$'],'fontsize',fontsize,'Interpreter','latex');
    ylabel('$PDF$','fontsize',fontsize,'Interpreter','latex');
    legend([h1 h2],{'Simu','Theo'},'fontsize',fontsize,'location','northeast','box','on');
    title(titlename1{isub},'fontsize',fontsize);
    H=get(gca,'position');
    set(gca,'position',[H(1)+0.02,H(2:end)]);

end
% subplot(2,2,1);
% %% data
% counts=histcounts(thetap_dat,countx)';
% ploty=counts/length(thetap_dat);
% plot(plotx,ploty,'k.-');hold on;
% %% theo
% counts=histcounts(thetap,countx)';
% ploty=counts/length(thetap);
% plot(plotx,ploty,'r.--');
% xlim([-90 90]);
% ylim([0 yylim]);
% set(gca,'xTick',-xxlim:30:xxlim);
% 
% plot([max(plinc) max(plinc)],[0 1],'k--');
% plot([-max(plinc) -max(plinc)],[0 1],'k--');
% xlabel('$\theta_{v_P}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% 
% subplot(2,2,2);
% counts=histcounts(phip_dat,countx_phi)';
% ploty=counts/length(phip_dat);
% plot(plotx_phi,ploty,'k.-');hold on;
% 
% counts=histcounts(phip,countx_phi)';
% ploty=counts/length(phip);
% plot(plotx_phi,ploty,'r.--');
% xlim([0 180]);
% ylim([0 yylim]);
% set(gca,'xTick',-xxlim_phi:30:xxlim_phi);
% xlabel('$\varphi_{v_P}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% 
% % plot([90+deltaang_max 90+deltaang_max],[0 1],'k--');
% % plot([90-deltaang_max 90-deltaang_max],[0 1],'k--');
% % 
% % plot([90+deltaang_min 90+deltaang_min],[0 1],'k--');
% % plot([90-deltaang_min 90-deltaang_min],[0 1],'k--');
% plot([90 90],[0 1],'k--');
% 
% subplot(2,2,3);
% counts=histcounts(thetav_dat,countx)';
% ploty=counts/length(thetav_dat);
% plot(plotx,ploty,'k.-');hold on;
% 
% counts=histcounts(thetav,countx)';
% ploty=counts/length(thetav);
% plot(plotx,ploty,'r.--');
% xlim([-xxlim xxlim]);
% ylim([0 yylim]);
% set(gca,'xTick',-xxlim:30:xxlim);
% xlabel('$\theta_{v_r}~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% 
% subplot(2,2,4);
% counts=histcounts(theta_dat,countx)';
% ploty=counts/length(theta_dat);
% plot(plotx,ploty,'k.-');hold on;
% 
% counts=histcounts(theta,countx)';
% ploty=counts/length(theta);
% plot(plotx,ploty,'r.--');
% xlim([-xxlim xxlim]);
% ylim([0 0.1]);
% set(gca,'xTick',-xxlim:30:xxlim);
% xlabel('$\theta_0~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% 
% 
