%% CE location distribution
clear;
ffname='RealPlutinosNpl';
fname={'1999CE119_1Gyr_40pl';'2001FU172_1Gyr_40pl';'1999CE119&2006RJ103_1Gyr_40pl';'2001FU172&2006RJ103_1Gyr_40pl'};
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

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
    %rectangle('position',[-100,-100,200,200],'facecolor','k');
    plot(Rt*cos(xt),Rt*sin(xt),'r-','linewidth',2);

    plot([-hillmax hillmax],[0 0],'r--');
    plot([0 0],[-hillmax hillmax],'r--');
    plotx=(PV_rlt_con(:,1).^2+PV_rlt_con(:,2).^2).^(1/2)/hill;
    ploty=PV_rlt_con(:,3)/hill;
    
    di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/di_record_inout.txt'));
    miu=1.94148e-12;
    absdi=abs(di_record_inout.*hnorm./Rnorm./coswf/miu);
    if isub==1
        maxdi=max(absdi);  
    end
    colorindex=abs(log(absdi/maxdi));
    colorindexNormal=(-erf(colorindex/3-2)+1)/2;
    data=zeros(length(plotx),3);
    data(:,1)=plotx;
    data(:,2)=ploty;
    data(:,3)=colorindexNormal;
    datasort=sortrows(data,3);
    color=[datasort(:,3) zeros(length(datasort),1) 1-datasort(:,3)];

    scatter(datasort(:,1),datasort(:,2),5,color,'filled');
    
    t=-pi/2:0.01:pi/2;
    r=0.01:0.01:3.5-0.001;
    [tt,rr]=meshgrid(t,r);
    zz=2*(1-(rr./hillmax).^2).^(1/2)./rr.*sin(tt);

    [xx,yy]=pol2cart(tt,rr);
    maxzz=max(max(zz));
    ev=exp(exp(0:0.1:log(log(1+maxzz)+1))-1)-1;
    ev=[ev -ev];
    ev=unique(ev);
    ev=sort(ev);
    contour(xx,yy,zz,ev,'w-.','linewidth',2);

    hold off;
    xlim([-hillmax hillmax]);
    ylim([-hillmax hillmax]);
    xlabel('$$\Vert(x_r,y_r)\Vert_2~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
    ylabel('$$z_r~(R_{H})$$','fontsize',fontsize,'Interpreter','latex');
    title(titlename{isub},'fontsize',fontsize);
    set(text(-3.0,1/9*hillmax,['$$N_{CE}= ',num2str(length(PV_record_tp)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','red');

end
