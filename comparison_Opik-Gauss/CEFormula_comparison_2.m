
clear;
aP=30;
at=40;et=0.25;it=30; % at=39.99999
% at=40.0;et=0.2500001;it=30;
% at=30.;et=0.00001;it=20;
OP=150;Ot=-30;wt=180;
mP=6e-9; % mP=3e-6;
ft_min=2.5e-5;
ft_max=0.25;
ft_max_trace=0.35; %deg
CEth=3.5;
Rth=CEth*aP.*(mP/3).^(1/3);
mPluto=6e-9;
Rmin=((mP/mPluto)^(1/3)*1188.3)/(Rth*1.5e8);

xlimSet=[4e-5,1];

% Pluto 半径是1188.3km 希尔半径是3.5*40*(1e-9/3)^(1/3)*1.5e8=1.5e7km
% 所以半径/希尔半径~=1e-4
% 用Pluto的平均密度来等比例缩放半径
% 所以如果1e-9(Pluto)质量的话，这里取f=1e-4就够了
% ft=[-exp(log(1e-5):0.08:log(ft_max))];
ft=-exp(log(ft_min):0.2:log(ft_max));
ft_trace=[-exp(log(ft_min):0.2:log(ft_max_trace)),exp(log(ft_min):0.2:log(ft_max_trace))];
% Begind=2;
% Endind=round(length(ft)/2)-10;
% at=39.9:0.001:40;

%% Constants


%% Detect independent variable
varL={'aP','OP','mP','at','et','it','Ot','wt','ft'};
fDetect=0;
for ivar=1:length(varL)
    var=varL{ivar};
    if length(eval(var))>1
        if fDetect
            error('Two or more arrays set!');
        end
        fDetect=1;
        fprintf('Array variable: %s\n', var);
        eval([var,'L=',var,';']);
        varSet=var;
    end
end
if ~fDetect
    error('No array set!');
end
clear var;

%% Initialization
len=length(eval(varSet));
% 0: Opik; 1: Gauss
di=zeros(len,1);di1=zeros(len,1);
de=zeros(len,1);de1=zeros(len,1);
da=zeros(len,1);da1=zeros(len,1);
Rr=zeros(len,1);
xbL=zeros(len,1);ybL=zeros(len,1);zbL=zeros(len,1);
x0L=zeros(len,1);y0L=zeros(len,1);z0L=zeros(len,1);
UxL=zeros(len,1);UyL=zeros(len,1);UzL=zeros(len,1);

%% Calculation
for ix=1:len
    eval([varSet,'=',varSet,'L(ix);']); 
    [di(ix),de(ix),da(ix),xb,yb,zb,sinPhi,cosPhi,x0,y0,z0,Ux,Uy,Uz] = ...,
        Fun_CEFormula_Opik(aP,OP,mP,at,et,it,Ot,wt,ft,CEth);
    x0L(ix)=x0*aP;y0L(ix)=y0*aP;z0L(ix)=z0*aP;
    UxL(ix)=Ux*aP;UyL(ix)=Uy*aP;UzL(ix)=Uz*aP;
    [di1(ix),de1(ix),da1(ix),Rr(ix)] = ...,
        Fun_CEFormula_Gauss(aP,OP,mP,at,et,it,Ot,wt,xb,yb,zb,sinPhi,cosPhi,CEth);
    xbL(ix)=xb*aP;ybL(ix)=yb*aP;zbL(ix)=zb*aP;
end

% %% Visualization-position
% %%% Visual interval
% figure;
% fontsize=20;
% set(gcf,'Position',[400,100,700,400],'color','w');
% 
% % lim=6e-2;
% % plot3(0,0,0,'b+','markersize',20);hold all;
% % axis([-lim lim -lim lim -lim lim]);axis square;
% % plot3(xbL(1),ybL(1),zbL(1),'r+');
% % plot3(xbL,ybL,zbL,'r.');grid on;
% % plot3(x0L(1),y0L(1),z0L(1),'k+');
% % plot3(x0L,y0L,z0L,'k.-');
% % % quiver3(x0L,y0L,z0L,UxL,UyL,UzL,'k')
% % quiver3(x0L,y0L,z0L,UxL/100,UyL/100,UzL/100,'k-',...,
% %     'autoscalefactor',0.9,'showarrowhead','off','maxheadsize',0.0001);
% % % quiver3(x0L,y0L,z0L,-UxL/10,-UyL/10,-UzL/10,'k-.',...,
% % %     'autoscalefactor',0.9,'showarrowhead','off','maxheadsize',0.0001);
% % hold off;
% 
% BottomRetainWidth0=0.1;
% LeftRetainWidth0=0.15;
% MiddleRetainWidth0=0.1;
% Height0=0.8;
% Width0=0.33;
% 
% len_trace=length(ft_trace);
% x0L_trace=zeros(len_trace,1);y0L_trace=zeros(len_trace,1);z0L_trace=zeros(len_trace,1);
% for ix2=1:len_trace
%     [x0L_trace(ix2),y0L_trace(ix2),z0L_trace(ix2)]=...,
%         Fun_CEFormula_EnterPos(aP,OP,at,et,it,Ot,wt,ft_trace(ix2));
% end
% x0L_trace=x0L_trace*aP;y0L_trace=y0L_trace*aP;z0L_trace=z0L_trace*aP;
% 
% subplot(1,2,1);
% xt=-2:0.1:2*pi;
% plot(Rth*cos(xt),Rth*sin(xt),'k--','linewidth',1);hold all;
% plot(0,0,'k.','markersize',50);
% quiver(y0L,x0L,UyL/100,UxL/100,'k-',...,
%     'autoscalefactor',0.9,'showarrowhead','on','maxheadsize',1);
% % plot(y0L(1),x0L(1),'k+');
% plot(y0L_trace,x0L_trace,'k-');
% plot(y0L,x0L,'k.','markersize',10);
% % plot(ybL(1),xbL(1),'r+');
% plot(ybL,xbL,'r.','markersize',10);
% % plot(0,0,'b+','markersize',20);
% set(gca,'position',[LeftRetainWidth0 BottomRetainWidth0 Width0 Height0]);
% pos0=get(gca,'position');
% pos0(1)=pos0(1)+pos0(3)/2;pos0(2)=pos0(2)+pos0(4)/2;
% pos0(3)=0;pos0(4)=-pos0(4)/4;
% annotation('arrow','position',pos0,'linestyle','--');
% boxHeight=0.05;
% annotation('textbox','position',[pos0(1)+pos0(3) pos0(2)+pos0(4)-boxHeight/2 0.1 boxHeight],...,
%     'string','Sun','linestyle','none','fontsize',fontsize,'Interpreter','latex');
% grid on;
% % axis([-lim lim -lim lim]);
% axis equal;axis square;
% xlabel('$y~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$x~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
% hold off;
% 
% subplot(1,2,2);
% plot(Rth*cos(xt),Rth*sin(xt),'k--','linewidth',1);hold all;
% plot(0,0,'k.','markersize',50);
% quiver(y0L,z0L,UyL/100,UzL/100,'k-',...,
%     'autoscalefactor',0.9,'showarrowhead','on','maxheadsize',1);
% % plot(y0L(1),z0L(1),'k+');
% plot(y0L_trace,z0L_trace,'k-');
% plot(y0L,z0L,'k.','markersize',10);
% % plot(ybL(1),zbL(1),'r+');
% plot(ybL,zbL,'r.','markersize',10);
% set(gca,'position',[LeftRetainWidth0+MiddleRetainWidth0+Width0 BottomRetainWidth0 Width0 Height0]);
% % plot(0,0,'b+','markersize',20);
% grid on;
% % axis([-lim lim -lim lim]);
% axis equal;axis square;
% xlabel('$y~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$z~\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');
% hold off;

%% Visualization-effect
figure;

set(gcf,'Position',[400,100,700,400],'color','w');
fontsize=20;
markersize=8;

% before abs, check sign
if ~all(di.*di1) || ~all(de.*de1) || ~all(da.*da1)
    error('Sign not consistent!');
end
% abs for log
% di=abs(di);di1=abs(di1);
% de=abs(de);de1=abs(de1);
% da=abs(da);da1=abs(da1);

%% ft as argument
% subplot(Nsub,1,1)
% loglog(ft,Rr,'k.');
% xlabel('$f\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$R_0/R_{th}$','fontsize',fontsize,'Interpreter','latex');
% 
% subplot(Nsub,1,2)
% loglog(ft,di,'k.');hold on;
% loglog(ft,di1,'r.');
% xlabel('$f\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\Delta i\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% 
% subplot(Nsub,1,3)
% loglog(ft,de,'k.');hold on;
% loglog(ft,de1,'r.');
% xlabel('$f\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\Delta e$','fontsize',fontsize,'Interpreter','latex');
% 
% subplot(Nsub,1,4)
% loglog(ft,da,'k.');hold on;
% loglog(ft,da1,'r.');
% xlabel('$f\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
% ylabel('$\Delta a\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');

%% Rr as argument
% 插值来把f写在上坐标轴把

BottomRetainWidth=0.15;
LeftRetainWidth=0.2;
Height=0.25;
Width=0.6;

dName={'di','de'}; %,'da'};
ylabelList={'$\Delta i~\mathrm{(DEG)}$','$\Delta e$'}; %,'$\Delta a~\mathrm{(AU)}$'};
% Nsub=length(dName)+1;
Nsub=length(dName);

%% Sub2-4
for iplot=1:Nsub
    
    axes('position',[LeftRetainWidth BottomRetainWidth+(Nsub-iplot)*Height Width Height]);

    dx=eval(dName{iplot});
    dx1=eval([dName{iplot},'1']);
    if strcmp(dName{iplot},'di')
        dx=-dx;
        dx1=-dx1;
    end
    rltErr=(abs(dx1)-abs(dx))./abs(dx);
    rltErrPos=rltErr;rltErrNeg=rltErr;
    lastNegInd=find(rltErrPos<0,1,'last');
    rltErrPos(lastNegInd)=-rltErrPos(lastNegInd);
    rltErrPos(rltErrPos<0)=nan;
    rltErrNeg(rltErrNeg>0)=nan;
    [ax,h1,h2]=plotyy(Rr,[dx,dx1],Rr,[rltErrPos,-rltErrNeg],@loglog,@loglog);hold all;
    %% axe properties
    set(ax,'xlim',xlimSet);
    set(ax(1),'ycolor','k');
    set(ax(2),'ycolor','k');
    set(get(ax(1),'Ylabel'),'String',ylabelList{iplot},'fontsize',fontsize,'Interpreter','latex');
    set(get(ax(2),'Ylabel'),'String','$|Rlt~Err|$','fontsize',fontsize,'Interpreter','latex');
    %%% The tick overlap problem
    linkaxes(ax,'x');
%     set(ax(1),'Box','on');
%     set(ax(2),'Box','off');
    set(ax(2), 'XTickLabel','','XAxisLocation','Top');
    %%% xlabel
    switch iplot
        case 1
            set(ax,'Box','on');
            set(ax,'xticklabel',[]);
            
        case 2
            set(ax,'Box','off');
            pos=get(ax(1),'position');
            annotation('line',[pos(1) pos(1)+pos(3)],[pos(2)+pos(4) pos(2)+pos(4)]);

            xlabel('$\gamma_R$','fontsize',fontsize,'Interpreter','latex');
    end
    %%% left 
    yylim=get(ax(1),'ylim');
    switch iplot
        case 1
            set(ax(1),'ytick',10.^(log10(yylim(1))+1:log10(yylim(2))));
        case 2
            set(ax(1),'ytick',10.^(log10(yylim(1)):log10(yylim(2))-1));
    end
    yylim=get(ax(2),'ylim');
    switch iplot
        case 1
            set(ax(2),'ytick',10.^(log10(yylim(1))+1:log10(yylim(2))));
        case 2
            set(ax(2),'ylim',[1e-9,1e0]);
            yylim=get(ax(2),'ylim');
            set(ax(2),'ytick',10.^(log10(yylim(1)):log10(yylim(2))-1));
            %         set(ax(2),'yminortick','on');
    end
    set(ax(1),'yticklabelmode','auto');
    set(ax(2),'yticklabelmode','auto');
    
    %%% neg yticklabel
    if strcmp(dName{iplot},'di')
        yticklabel=get(ax(1),'yticklabel');
        for ic=1:length(yticklabel)
            yticklabel{ic}=['-',yticklabel{ic}];
        end
        set(ax(1),'yticklabel',yticklabel);
    end
%     yticklabel2=get(ax(2),'yticklabel');
%     for ic=1:length(yticklabel2)
%         yticklabel2{ic}=['-',yticklabel2{ic}];
%     end
%     set(ax(2),'yticklabel',yticklabel2);
%     box on;
    %% line properties
    set(h1(1),'linestyle','-','color','b','linewidth',2.0);
    set(h1(2),'linestyle','none','marker','.','color','r','markersize',12);
    set(h2(1),'linestyle','--','color',[1 0.5 0],'linewidth',1.5);
    set(h2(2),'linestyle','--','color',[0 1 0],'linewidth',1.5);
    if strcmp(dName{iplot},'di')
        legend([h1(1),h1(2)],{'Opik','Gauss'},'fontsize',15,'Interpreter','latex','location','southwest');
    elseif strcmp(dName{iplot},'de')
        legend([h2(1),h2(2)],{'$+$','$-$'},'fontsize',15,'Interpreter','latex','location','southwest');
    end
    %% annotation
    %%% interior shadow
    xxlim=get(gca,'xlim');yylim=get(gca,'ylim');
    patch([xxlim(1) xxlim(1) Rmin Rmin],...,
        [yylim(1) yylim(2) yylim(2) yylim(1)],...,
        'k','facealpha',0.70,'edgecolor','none');
    plot([Rmin Rmin],yylim,'k-');
    %%% 0.01 line
    xxlim=get(ax(2),'xlim');yylim=get(ax(2),'ylim');
    line(xxlim,[0.01 0.01],'parent',ax(2),'color','k','linestyle','--');
    pos=get(ax(2),'position');
    arrowPos=[pos(1)+pos(3)*0.99,...,
              pos(2)+pos(4)/(log(yylim(2))-log(yylim(1)))*(log(0.01)-log(yylim(1))),...,
              pos(3)*0.01 0];
    annotation('arrow','position',arrowPos,'linestyle','-','color','k');
    hold off;

end

