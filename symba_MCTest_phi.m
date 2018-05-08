if ~exist('phi1_1','var')
  
clear;

ffname1='symbaRealPlutinosNpl_fast';
ffname2='symbaRanPlutinosNpl';
fname1={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
titlename1={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};
Dir='ServerMount';

for isub=1:length(fname1)
    
    fname=fname1{isub};
    titlename=titlename1{isub};
    disp(titlename);

    for iplot=1:2
        ffname=eval(['ffname',num2str(iplot)]);
        disp(ffname);

        r2hill_record=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt']);
        r2hill_mean=mean(r2hill_record);
        hill=sqrt(r2hill_mean);
        PV_record_pl=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt']);
        PV_record_tp=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt']);
        AE_record_pl=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_pl.txt']);
        AE_record_tp=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/AE_record_tp.txt']);
        
        if iplot==1
            ierror=load(['~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt']);
            PV_record_pl(ierror,:)=[];
            PV_record_tp(ierror,:)=[];
            AE_record_pl(ierror,:)=[];
            AE_record_tp(ierror,:)=[];
        end
        
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
              
        x=P_rlt_con(:,1);
        y=P_rlt_con(:,2);
        
        phi=circ_rad2ang(angle(x+y*1i));
        phi=phi+180;
        eval(['phi',num2str(iplot),'_',num2str(isub),'=phi;']);
        
    end
end

end

%% Plot
fontsize=15;
N=50;
yylim=0.02;
xxlim=[0 360];
colorList={'k','r'};

phix=0:N;
phix=phix';
dphi=(xxlim(2)-xxlim(1))/N;
phix=xxlim(1)+phix*dphi;
plotx=(phix(2:end)-dphi/2);

figure;
set(gcf,'Position',[400,100,700,500],'color','w');
for isub=1:length(fname1)
    
    subplot(2,2,isub);
    plot(0,0,'w');hold all;

    for iplot=1:2
                            
        eval(['phi=phi',num2str(iplot),'_',num2str(isub),';']);
        color=colorList{iplot};
        
        countx=histcounts(phi,phix)'/dphi;
        ploty=countx/length(phi);
        
        %bar(plotx,ploty,'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        h=plot(plotx,ploty,[color,'.-'],'linewidth',1.3,'markersize',10);
        eval(['h',num2str(iplot),'=h;']);
        
        xlim(xxlim);
        ylim([0 yylim]);
        
        if iplot==1
            plot([90 90],[0 yylim],[color,'--']);
            plot([270 270],[0 yylim],[color,'--']);
        end
    end
    hold off;
    
    legend([h1 h2],{'Num','MC'},'box','on','fontsize',fontsize);
    title(titlename,'fontsize',fontsize);

    set(gca,'xTick',xxlim(1):30:xxlim(2));
    xlabel('$\varphi_0~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$PDF$','fontsize',fontsize,'Interpreter','latex');
    H=get(gca,'position');
    set(gca,'position',[H(1)+0.02,H(2:end)]);

end
