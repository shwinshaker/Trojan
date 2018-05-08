
clear;
aP=30;at=40;et=0.25;it=30;
OP=150;Ot=-30;wt=180;
mPL=exp(log(1e-9):0.2:log(1e-5));
Rr=1e-2;

xlimSet=[5e-10,2e-5];

%% Constants
CEth=3.5;
Rth=CEth*aP*(mPL/3).^(1/3);
mPluto=6e-9;mEarth=3e-6;mMercury=0.055*mEarth;
Rmin=((mPL/mPluto).^(1/3)*1188.3)./(Rth*1.5e8);

%% Rr-ft fit
%%% Rb和ft之间的关系不受mP影响(计算Rb无需mP)
%%% 所以在此做一次拟合就行，随后其他的mP就按立方根缩放就行了
%%% st->Standard
mP_st=6e-9;% Rth_st=3.5*(mP_st/3).^(1/3);
ft_max_st=2.5e-2;ft_min_st=2.5e-5;
ft_st=-exp(log(ft_min_st):0.2:log(ft_max_st))';
RbL=zeros(length(ft_st),1);
for ix=1:length(ft_st)
    [~,~,~,xb,yb,zb,~,~,~,~,~,~,~,~] = ...,
        Fun_CEFormula_Opik(aP,OP,mP_st,at,et,it,Ot,wt,ft_st(ix),CEth);
    RbL(ix)=sqrt(xb^2+yb^2+zb^2);
end
P=polyfit(RbL,ft_st,1);
fprintf('\nk=%.7f b=%.7f\n',P(1),P(2));

%% determine respective f to keep Rr constant for all mP
ftL=polyval(P,Rr*Rth/aP);
[~,~,~,xb,yb,zb,~,~,~,~,~,~,~,~] = ...,
    Fun_CEFormula_Opik(aP,OP,mPL(1),at,et,it,Ot,wt,ftL(1),CEth);
Rr_1=sqrt(xb^2+yb^2+zb^2)/(Rth(1)/aP);
fprintf('Rr_1:%.7f\n',Rr_1);

%% Initialization
len=length(mPL);
% 0: Opik; 1: Gauss
di=zeros(len,1);di1=zeros(len,1);
de=zeros(len,1);de1=zeros(len,1);
da=zeros(len,1);da1=zeros(len,1);
xbL=zeros(len,1);ybL=zeros(len,1);zbL=zeros(len,1);
x0L=zeros(len,1);y0L=zeros(len,1);z0L=zeros(len,1);
UxL=zeros(len,1);UyL=zeros(len,1);UzL=zeros(len,1);

%% Calculation
for ix=1:len
    mP=mPL(ix);
    ft=ftL(ix);
    [di(ix),de(ix),da(ix),xb,yb,zb,sinPhi,cosPhi,x0,y0,z0,Ux,Uy,Uz] = ...,
        Fun_CEFormula_Opik(aP,OP,mP,at,et,it,Ot,wt,ft,CEth);
    x0L(ix)=x0*aP;y0L(ix)=y0*aP;z0L(ix)=z0*aP;
    UxL(ix)=Ux*aP;UyL(ix)=Uy*aP;UzL(ix)=Uz*aP;
    [di1(ix),de1(ix),da1(ix),Rr1] = ...,
        Fun_CEFormula_Gauss(aP,OP,mP,at,et,it,Ot,wt,xb,yb,zb,sinPhi,cosPhi,CEth);
    if abs(Rr1-Rr)>0.0001
        fprintf('Rr1:%.7f\n',Rr1);
        error('Bad fitting!');
    end
    xbL(ix)=xb*aP;ybL(ix)=yb*aP;zbL(ix)=zb*aP;
end

%% Visualization-effect
figure;
set(gcf,'Position',[400,100,700,400],'color','w');
fontsize=20;
markersize=8;

% before abs, check sign
if ~all(di.*di1) || ~all(de.*de1) || ~all(da.*da1)
    error('Sign not consistent!');
end

%% Rr as argument
% 插值来把f写在上坐标轴把
BottomRetainWidth=0.15;
LeftRetainWidth=0.2;
Height=0.25;
Width=0.6;

dName={'di','de'}; %,'da'};
Nsub=length(dName);
ylabelList={'$\Delta i~\mathrm{(DEG)}$','$\Delta e$'}; % ,'$\Delta a~\mathrm{(AU)}$'};

Nant=1;

%% Sub1,2
for iplot=1:2
    
%     subplot(Nsub,1,iplot)
    switch iplot
        case 1
            axes('position',[LeftRetainWidth BottomRetainWidth+(Nsub-iplot)*Height Width Height/4*(4+Nant+1)]);

        case 2 
            axes('position',[LeftRetainWidth BottomRetainWidth+(Nsub-iplot)*Height Width Height]);
    end
    % set(ax,'position',[LeftRetainWidth BottomRetainWidth+(Nsub-iplot-1)*Height Width Height]);
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
    [ax,h1,h2]=plotyy(mPL,[dx,dx1],mPL,[rltErrPos,-rltErrNeg],@loglog,@loglog);hold all;
    %% axe properties
    set(ax,'xlim',xlimSet);
    set(get(ax(1),'Ylabel'),'String',ylabelList{iplot},'fontsize',fontsize,'Interpreter','latex');
    set(get(ax(2),'Ylabel'),'String','$|Rlt~Err|$','fontsize',fontsize,'Interpreter','latex');
    set(ax(1),'ycolor','k');
    set(ax(2),'ycolor','k');
    %%% The tick overlap problem
    linkaxes(ax,'x');
    % set(ax(1),'Box','off');
    % set(ax(2),'Box','off');
    set(ax(2), 'XTickLabel','','XAxisLocation','Top');
    pos=get(ax(1),'position');
    % annotation('line',[pos(1) pos(1)+pos(3)],[pos(2)+pos(4) pos(2)+pos(4)]);
    %%% xlabel
    if iplot==2
        xlabel('$m_P/m_\odot$','fontsize',fontsize,'Interpreter','latex');
    else
        set(ax,'xticklabel',[]);
        title(['$\gamma_R=',sprintf('%.2f',Rr),'$'],'fontsize',fontsize,'Interpreter','latex');
    end
    %% left plot
    switch iplot
        case 1;
            yylim=get(ax(1),'ylim');
            set(ax(1),'ylim',[yylim(1) 10^(log10(yylim(2))+Nant)]);
            yylim=get(ax(1),'ylim');
            set(ax(1),'ytick',10.^(floor(log10(yylim(1)))+1:ceil(log10(yylim(2)))));
        case 2;
            set(ax(1),'ylim',[1e-4 1]);
            yylim=get(ax(1),'ylim');
            set(ax(1),'ytick',10.^(floor(log10(yylim(1))):ceil(log10(yylim(2)))-1));
    end
    %% right plot
    switch iplot
        case 1;
            set(ax(1), 'XTickLabel','','XAxisLocation','Top');
            yylim=get(ax(2),'ylim');
            set(ax(2),'ylim',[yylim(1) 10^(log10(yylim(2))+Nant+1)]);
            yylim=get(ax(2),'ylim');
            set(ax(2),'ytick',10.^(floor(log10(yylim(1)))+1:ceil(log10(yylim(2)))));
        case 2;
            set(ax(2),'ylim',[1e-6 1e-2]);
            yylim=get(ax(2),'ylim');
            set(ax(2),'ytick',10.^(floor(log10(yylim(1))):ceil(log10(yylim(2)))-1));
    end
    set(ax(1),'yticklabelmode','auto');
    set(ax(2),'yticklabelmode','auto');
    %% negative label
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
    %% line properties
    set(h1(1),'linestyle','-','color','b','linewidth',2);
    set(h1(2),'linestyle','none','marker','.','color','r','markersize',12);
    set(h2(1),'linestyle','-','linewidth',1.5,'color',[1 0.5 0]);
    set(h2(2),'linestyle','-','linewidth',1.5,'color',[0 1 0]);
    if strcmp(dName{iplot},'di')
        legend([h1(1),h1(2)],{'Opik','Gauss'},'fontsize',15,'Interpreter','latex','location','southeast');
    elseif strcmp(dName{iplot},'de')
        legend([h2(1),h2(2)],{'$+$','$-$'},'fontsize',15,'Interpreter','latex','location','southeast');
    end
    %% annotation
    % patch([xxlim(1) xxlim(1) Rmin Rmin],...,
    %     [yylim(1) yylim(2) yylim(2) yylim(1)],...,
    %     'k','facealpha',0.70,'edgecolor','none');
    % plot([Rmin Rmin],yylim,'k-');
    xxlim=get(ax(2),'xlim');yylim=get(ax(2),'ylim');
    if yylim(1)<=0.01 && yylim(2)>=0.01
        line(xxlim,[0.01 0.01],'parent',ax(2),'color','k','linestyle','--');
        pos=get(ax(2),'position');
        arrowPos=[pos(1)+pos(3)*0.99,...,
                  pos(2)+pos(4)/(log(yylim(2))-log(yylim(1)))*(log(0.01)-log(yylim(1))),...,
                  pos(3)*0.01 0];
        annotation('arrow','position',arrowPos,'linestyle','-','color','k');
    end
    if iplot==1
        arrowVerticalLengthRatio=0.12;
        pos=get(gca,'position');
        xxlim=get(gca,'xlim');yylim=get(gca,'ylim');
        %%% Pluto
        arrowPos=[pos(1)+pos(3)/(log(xxlim(2))-log(xxlim(1)))*(log(mPluto)-log(xxlim(1))) ...,
                  pos(2)+pos(4)*(1-arrowVerticalLengthRatio) ...,
                  0 ...,
                  pos(4)*arrowVerticalLengthRatio];
        annotation('arrow','position',arrowPos,'linestyle','-','color','red');
        annotation('textbox','position',[arrowPos(1) arrowPos(2) -0.03 -0.05],...,
            'string','Pluto','linestyle','none','fontsize',fontsize,'Interpreter','latex','color','r');
        %%% Earth
        arrowPos=[pos(1)+pos(3)/(log(xxlim(2))-log(xxlim(1)))*(log(mEarth)-log(xxlim(1))) ...,
                  pos(2)+pos(4)*(1-arrowVerticalLengthRatio) ...,
                  0 ...,
                  pos(4)*arrowVerticalLengthRatio];
        annotation('arrow','position',arrowPos,'linestyle','-','color','red');
        annotation('textbox','position',[arrowPos(1) arrowPos(2) -0.03 -0.05],...,
            'string','Earth','linestyle','none','fontsize',fontsize,'Interpreter','latex','color','r');
        %%% Mercury
        arrowPos=[pos(1)+pos(3)/(log(xxlim(2))-log(xxlim(1)))*(log(mMercury)-log(xxlim(1))) ...,
                  pos(2)+pos(4)*(1-arrowVerticalLengthRatio) ...,
                  0 ...,
                  pos(4)*arrowVerticalLengthRatio];
        annotation('arrow','position',arrowPos,'linestyle','-','color','red');
        annotation('textbox','position',[arrowPos(1) arrowPos(2) -0.03 -0.05],...,
            'string','Mercury','linestyle','none','fontsize',fontsize,'Interpreter','latex','color','r');
    end
    hold off;

end


