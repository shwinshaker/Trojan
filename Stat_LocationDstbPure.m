%% CE location distribution
clear;
% ffname='RealPlutinosNpl';
% fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

ffname='RanPlutinosNpl';
fname={'1999CE119_1Gyr_40pl_ran';'2001FU172_1Gyr_40pl_ran';'1999CE119&2006RJ103_1Gyr_40pl_ran';'2001FU172&2006RJ103_1Gyr_40pl_ran'};

hillmax=3.5;
fontsize=15;

figure(1);
set(gcf,'Position',[400,100,700,700],'color','w');


for isub=1:4
    subplot(2,2,isub);
    
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_tp.txt'));
    AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_pl.txt'));
    AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_tp.txt'));
    r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/r2hill_record.txt'));
    r2hill_mean=mean(r2hill_record);
    hill=sqrt(r2hill_mean);

    PV_rlt=PV_record_pl-PV_record_tp;
    PV_rlt=PV_rlt(:,1:3);
    dis=(PV_rlt(:,1).^2+PV_rlt(:,2).^2+PV_rlt(:,3).^2).^(1/2);
    Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
    hnorm=((2*pi)^2/(365.25)^2*AE_record_tp(:,1).*(1-AE_record_tp(:,2).^2)).^(1/2);
    sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
    coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));


    PV_rlt_con=zeros(length(PV_rlt),3);
    PV_rlt_con_norm=zeros(length(PV_rlt),1);
    PV_tp_norm=zeros(length(PV_rlt),1);
    
    %% verify
    Ptp=PV_record_tp(:,1:3);
%     Vtp=PV_record_tp(:,4:6);
%     h=[Ptp(:,2).*Vtp(:,3)-Ptp(:,3).*Vtp(:,2),Ptp(:,3).*Vtp(:,1)-Ptp(:,1).*Vtp(:,3),Ptp(:,1).*Vtp(:,2)-Ptp(:,2).*Vtp(:,1)]; 
    Ptp_con=zeros(length(PV_rlt),3);
    
    for i=1:length(PV_rlt)
     P1=[cosd(AE_record_tp(i,4)) sind(AE_record_tp(i,4)) 0; ...,
        -sind(AE_record_tp(i,4)) cosd(AE_record_tp(i,4)) 0; ...,
        0 0 1];

     P2=[1 0 0;...,
        0 cosd(AE_record_tp(i,3)) sind(AE_record_tp(i,3)); ...,
        0 -sind(AE_record_tp(i,3)) cosd(AE_record_tp(i,3))];

     P3=[coswf(i) sinwf(i) 0;...,
        -sinwf(i) coswf(i) 0;...,
         0 0 1];
     
     %% verify
     Ptp_con(i,:)=(P3*P2*P1*Ptp(i,:)')';
     PV_rlt_con(i,:)=(P3*P2*P1*PV_rlt(i,:)')';
     PV_rlt_con_norm(i)=norm(PV_rlt_con(i,:));
     PV_tp_norm(i)=norm(PV_record_tp(i,:));
    end
    
    %
%     I=zeros(234, 123);
%     imshow(I);hold all;
    
    xt=-pi:0.1:2*pi;
    Rt=hillmax;
    
    plot(Rt*cos(xt),Rt*sin(xt),'w-','linewidth',2);axis square;hold all;
    plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);

    plot([-hillmax hillmax],[0 0],'r--');
    plot([0 0],[-hillmax hillmax],'r--');
    
    %% x z
%     plotx=(PV_rlt_con(:,1).^2+PV_rlt_con(:,2).^2).^(1/2)/hill;
%     ploty=PV_rlt_con(:,3)/hill;
    
    %% x y
    plotx=PV_rlt_con(:,1)/hill;
    ploty=PV_rlt_con(:,2)/hill;
    
    plot(plotx,ploty,'k.','markersize',3);

    hold off;
    xlim([-hillmax hillmax]);
    ylim([-hillmax hillmax]);
    xlabel('$$x_r~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$$y_r~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
    title(titlename{isub},'fontsize',fontsize);
%     set(text(-3.0,1/9*hillmax,['$$N_{CE}= ',num2str(length(PV_record_tp)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');

end
