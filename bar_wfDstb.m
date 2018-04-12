clear;
clear;
ffname='RealPlutinos';
fname={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

fontsize=15;
N=50;
yylim=0.1;

set(gcf,'Position',[400,100,700,700],'color','w');

for isub=1:4
    PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_pl.txt'));
    PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/PV_record_tp.txt'));
    AE_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_pl.txt'));
    AE_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/AE_record_tp.txt'));
    
    Rnorm=(PV_record_tp(:,1).^2+PV_record_tp(:,2).^2+PV_record_tp(:,3).^2).^(1/2);
    %hnorm=((2*pi)^2/(365.25)^2*AE_record_tp(:,1).*(1-AE_record_tp(:,2).^2)).^(1/2);
    %h=cross(PV_record_tp(:,1:3),PV_record_tp(:,4:6));
    %hnorm2=(h(:,1).^2+h(:,2).^2+h(:,3).^2).^(1/2);
    %Rdot=
    sinwf=PV_record_tp(:,3)./(Rnorm.*sind(AE_record_tp(:,3)));
    coswf=secd(AE_record_tp(:,4)).*(PV_record_tp(:,1)./Rnorm+sind(AE_record_tp(:,4)).*sinwf.*cosd(AE_record_tp(:,3)));
    wf=mod(angle(coswf+sinwf*1i)/pi*180,360);
    %ae=AE_record_tp(:,1).*(1-AE_record_tp(:,2).^2);
    cosw=cosd(AE_record_tp(:,5));
    sinw=sind(AE_record_tp(:,5));
    cosf=coswf.*cosw+sinwf.*sinw;
    sinf=sinwf.*cosw-coswf.*sinw;
    f=mod(angle(cosf+sinf*1i)/pi*180,360);
    
    Rnormpl=(PV_record_pl(:,1).^2+PV_record_pl(:,2).^2+PV_record_pl(:,3).^2).^(1/2);
    sinwfpl=PV_record_pl(:,3)./(Rnormpl.*sind(AE_record_pl(:,3)));
    coswfpl=secd(AE_record_pl(:,4)).*(PV_record_pl(:,1)./Rnormpl+sind(AE_record_pl(:,4)).*sinwfpl.*cosd(AE_record_pl(:,3)));
    wfpl=mod(angle(coswfpl+sinwfpl*1i)/pi*180,360);
    coswpl=cosd(AE_record_pl(:,5));
    sinwpl=sind(AE_record_pl(:,5));
    cosfpl=coswfpl.*coswpl+sinwfpl.*sinwpl;
    sinfpl=sinwfpl.*coswpl-coswfpl.*sinwpl;
    fpl=mod(angle(cosfpl+sinfpl*1i)/pi*180,360);
    
    subplot(2,2,isub);
%     plot(sind(AE_record_tp(:,5)),sinwf,'b.','markersize',0.5);hold on;
%     axis square;
%     plot(sind(AE_record_pl(:,5)),sinwfpl,'r.','markersize',0.5);
%     plot(AE_record_tp(:,2),sinwf,'b.','markersize',0.5);hold on;
%     axis square;
%     plot(AE_record_pl(:,2),sinwfpl,'r.','markersize',0.5);
    plot(f,wf,'r.','markersize',0.5);hold on;
    axis square;
    plot(fpl,wfpl,'k.','markersize',0.5);
%     plot(Rnorm,wf,'b.','markersize',0.5);hold on;
%     axis square;
%     plot(Rnormpl,wfpl,'r.','markersize',0.5);
%     plot(AE_record_tp(:,6),wf,'b.','markersize',0.5);hold on;
%     axis square;
%     plot(AE_record_pl(:,6),wfpl,'r.','markersize',0.5);



    xlim([0 360]);
    ylim([0 360]);
    set(gca,'xTick',0:30:360);
    set(gca,'yTick',0:30:360);


%     %%Mean
%     MeanAzm=mean(theta)/pi*180;
%     VarAzm=var(theta);
%     plot([MeanAzm MeanAzm],[0 yylim],'k--');

    hold off;
%     %set(text(20,8/9*yylim,['$$N_{CE}= ',num2str(length(theta)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     set(text(-20,8/9*yylim,['$$\overline{\theta}= ',num2str(MeanAzm,4),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');
%     title(titlename{isub},'fontsize',fontsize);
% 
%     xlabel('$\theta~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
%     ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

end




