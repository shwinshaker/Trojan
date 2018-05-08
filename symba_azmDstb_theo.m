
if ~exist('theta_dat1','var')
    
    clear;
    %% data
    ffname='symbaRealPlutinosNpl_fast';
    %% ffname='symba_RealPlutinos';
    ffname_theo='symba_RealPlutinos';
    fnameL={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
    titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
    % suf='40pl';
    
    inda=3;
    
    for ipl=1:4
        % fname=[fnameL{ipl},'_1Gyr_',suf];
        fname=fnameL{ipl};
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
        clear i;
        
        if isreal(i)
            error('Fatal error! i is assigned!');
        end
        
        %% Theta_P
        xv=(Vp_con(:,1).^2+Vp_con(:,2).^2).^(1/2);
        yv=Vp_con(:,3);
        thetap_dat=circ_rad2ang(angle(xv+yv*i));
        eval(['thetap_dat',num2str(ipl),'=thetap_dat;']);
        
        %% Phi_P
        phip_dat=circ_rad2ang(angle(Vp_con(:,1)+Vp_con(:,2)*i));
        eval(['phip_dat',num2str(ipl),'=phip_dat;']);
        
        %% phi_vPvT
        phivpvt_dat=acosd(sum(Vp_con.*Vt_con,2)./(sum(Vp_con.*Vp_con,2).^(1/2).*sum(Vt_con.*Vt_con,2).^(1/2)));
        eval(['phivpvt_dat',num2str(ipl),'=phivpvt_dat;']);
        
        %% Theta_v
        xv=(V_rlt_con(:,1).^2+V_rlt_con(:,2).^2).^(1/2);
        yv=V_rlt_con(:,3);
        thetav_dat=circ_rad2ang(angle(xv+yv*i));
        eval(['thetav_dat',num2str(ipl),'=thetav_dat;']);
        
        %% Phi_v
        xv=V_rlt_con(:,1);
        yv=V_rlt_con(:,2);
        phiv_dat=circ_rad2ang(angle(xv+yv*i));
%         phiv_dat(xv<0 & yv<0)=phiv_dat(xv<0 & yv<0)+180;
%         phiv_dat(xv<0 & yv>0)=phiv_dat(xv<0 & yv>0)+180;
        eval(['phiv_dat',num2str(ipl),'=phiv_dat;']);
        
        %% Theta
        xv=(P_rlt_con(:,1).^2+P_rlt_con(:,2).^2).^(1/2);
        yv=P_rlt_con(:,3);
        theta_dat=circ_rad2ang(angle(xv+yv*i));
        eval(['theta_dat',num2str(ipl),'=theta_dat;']);
        
        %% Phi
        xv=P_rlt_con(:,1);
        yv=P_rlt_con(:,2);
        phi_dat=circ_rad2ang(angle(xv+yv*i));
        phi_dat=phi_dat+180;
        eval(['phi_dat',num2str(ipl),'=phi_dat;']);
        
        %% lambda
        coswf(coswf>1)=1;
        coswf(coswf<-1)=-1;
        lambda_dat=circ_rad2ang(angle(coswf+sinwf*i));
        eval(['lambda_dat',num2str(ipl),'=lambda_dat;']);
        
        
        %% ki
        zeta=phiv_dat-phi_dat-90;
        tanki=tand(zeta)./sind(thetav_dat);
        sinki=sind(theta_dat)./cosd(thetav_dat);
        coski=sinki./tanki;
        ki_dat=circ_rad2ang(angle(coski+sinki*i));
        ki_dat=ki_dat+180;
        eval(['ki_dat',num2str(ipl),'=ki_dat;']);
        
        %% f
        cosf(cosf>1)=1;
        cosf(cosf<-1)=-1;
        f_dat=circ_rad2ang(angle(cosf+sinf*i));
        eval(['f_dat',num2str(ipl),'=f_dat;']);
        
    end

end
    
for ipl=1:4
    if ~exist(['theta',num2str(ipl)],'var')
        
        %% Theo theta Distribution
        disp([fnameL{ipl},'_theo']);
        %% get plel
        plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname_theo,'/',fnameL{ipl},'/plel.txt']);
        tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname_theo,'/',fnameL{ipl},'/tpel.txt']);
      
        pla=plel(:,inda);plamean=mean(pla);
        ple=plel(:,inda+1);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
        plinc=plel(:,inda+2);plincmean=mean(plinc);%plincstd=std(plinc);
        
        tpa=tpel(:,inda);tpamean=mean(tpa);
        tpe=tpel(:,inda+1);tpemean=mean(tpe);
        tpinc=tpel(:,inda+2);tpincmean=mean(tpinc);
        
        [theta,phi,thetap,phip,thetav,phiv,phivpvt,lambda,f,ki,~]=...,
            Fun_azmDstb_theo(plamean,plemean,plincmean,tpamean,tpincmean);
        
        eval(['thetap',num2str(ipl),'=thetap;']);
        eval(['phip',num2str(ipl),'=phip;']);
        eval(['phivpvt',num2str(ipl),'=phivpvt;']);
        eval(['thetav',num2str(ipl),'=thetav;']);
        eval(['phiv',num2str(ipl),'=phiv;']);
        eval(['theta',num2str(ipl),'=theta;']);
        eval(['phi',num2str(ipl),'=phi;']);
        eval(['lambda',num2str(ipl),'=lambda;']);
        eval(['f',num2str(ipl),'=f;']);
        eval(['ki',num2str(ipl),'=ki;']);
        
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
        xxlim1=-180;xxlim2=180;
    case 'theta'
        yylim=0.04;
        label='\theta_0';
    case 'phi'
        label='\varphi_0';
        yylim=0.02;
        xxlim1=0;xxlim2=360;
    case 'f'
        label='f_T';
        yylim=0.02;
        xxlim1=-180;xxlim2=180;
    case 'ki'
        label='ki';
        yylim=0.02;
        xxlim1=0;xxlim2=360;
    case 'lambda'
        xxlim1=-180;xxlim2=180;
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
    
    %% kuiper test
%     [h,k,K]=circ_kuipertest(circ_ang2rad(theta(theta>0)),abs(circ_ang2rad(theta(theta<0))));
%     fprintf('Kuiper test result: %i\n', h);
%     fprintf('Kuiper test result: %.4f', k);
%     fprintf('Kuiper test result: %.4f', K);

    % circ mean and confidence
    [cmean,cul,cll]=circ_mean(circ_ang2rad(var_dat));
    cmean=circ_rad2ang(cmean);
    cul=circ_rad2ang(cul);
    cll=circ_rad2ang(cll);
    conf=circ_confmean(circ_ang2rad(var_dat));
    conf=circ_rad2ang(conf);
    
%     plot([cll cll],[0 yylim],'k--');
%     plot([cul cul],[0 yylim],'k--');

    %% circ skewness
    skew=circ_skewness(circ_ang2rad(var_dat));
    
    set(text(xxlim1+0.05*(xxlim2-xxlim1),8/9*yylim,['$$\overline{\theta_0}=\left(',sprintf('%.4f',cmean),'\pm',sprintf('%.4f',conf),'\right)^\circ$$']),'Interpreter','latex','fontsize',fontsize/5*4,'color','red');
    set(text(xxlim1+0.05*(xxlim2-xxlim1),7/9*yylim,['$$Skew= ',sprintf('%.4f',skew),'$$']),'Interpreter','latex','fontsize',fontsize/5*4,'color','red');
    
%     h2=plot([xxlim1 xxlim2],[1/360 1/360],'r-','linewidth',1.3);

%     plot([90 90],[0 yylim],'k--');
%     plot([270 270],[0 yylim],'k--');
    plot([0 0],[0 yylim],'k--');
%     plot([180 180],[0 yylim],'k--');
    
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