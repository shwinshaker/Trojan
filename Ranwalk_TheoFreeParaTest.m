%% Theo theta Distribution

if ~(exist('aPL','var') && exist('aTL','var'))
    
    fnameList={'1999CE119','2001FU172'};
    Npl=length(fnameList);
    
    aPL=zeros(Npl);
    aTL=zeros(Npl);
    ePL=zeros(Npl);
    eTL=zeros(Npl);
    IPL=zeros(Npl);
    ITL=zeros(Npl);

    NCE_data=zeros(Npl);
    StdI_data=zeros(Npl);
    Stde_data=zeros(Npl);    
    AI_data=zeros(Npl);
    Ae_data=zeros(Npl);
        
    for ipl=1:Npl
        fname=fnameList{ipl};
    
        %% get plel
        plel=load(['~/Documents/ServerMount/LAB/CE_realp/RealPlutinos/',fname,'_1Gyr/plel.txt']);
        tpel=load(['~/Documents/ServerMount/LAB/CE_realp/RealPlutinos/',fname,'_1Gyr/tpel.txt']);
        
        aPL(ipl)=mean(plel(:,2));
        aTL(ipl)=mean(tpel(:,2));
        ePL(ipl)=mean(plel(:,3));
        eTL(ipl)=mean(tpel(:,3));
        IPL(ipl)=mean(plel(:,4))/180*pi;
        ITL(ipl)=mean(tpel(:,4))/180*pi;
        
        di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/RealPlutinos/',fname,'_1Gyr/di_record_inout.txt']);
        de_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/RealPlutinos/',fname,'_1Gyr/de_record_inout.txt']);
        di_record_inout=di_record_inout/180*pi;
        StdI_data(ipl)=std(di_record_inout);
        Stde_data(ipl)=std(de_record_inout);
        
        NCE_data(ipl)=length(di_record_inout);
        
        AI_data(ipl)=abs(sum(di_record_inout));
        Ae_data(ipl)=abs(sum(de_record_inout));
        %     AI_data=sqrt(2/pi*StdI_data^2*NCE_data);
        %     Ae_data=sqrt(2/pi*Stde_data^2*NCE_data);
    end
end

if ~(exist('NCE','var') && exist('MI','var') && exist('Me','var'))
    
    %% Theo
    NCE=zeros(Npl);
    for ipl=1:Npl
        NCE(ipl)=Fun_NCE(aPL(ipl),aTL(ipl),ePL(ipl),eTL(ipl),IPL(ipl),ITL(ipl),1);
    end
    
    stdList=0:0.2:5;
    Nls=length(stdList);
    
    AI=zeros(Nls,Npl);
    StdI=zeros(Nls,Npl);
    for ipl=1:Npl
        for il=1:Nls
            disp(stdList(il));
            [dinc,~,~,~,~,~,~,~,~,~,~]=Fun_diDstb_theo(100000,1,...,
                aPL(ipl),ePL(ipl),IPL(ipl)/pi*180,...,
                aTL(ipl),eTL(ipl),ITL(ipl)/pi*180,stdList(il));
            dinc=dinc/180*pi;
            StdI(il,ipl)=std(dinc);
            MI=StdI(il,ipl)^2*NCE(ipl);
            AI(il,ipl)=sqrt(2/pi*MI);
        end
    end

    Ae=zeros(Nls,Npl);
    Stde=zeros(Nls,Npl);
    for ipl=1:Npl
        for il=1:Nls
            disp(stdList(il));
            [~,de,~,~,~,~,~,~,~,~,~]=Fun_diDstb_theo(100000,1,...,
                aPL(ipl),ePL(ipl),IPL(ipl)/pi*180,...,
                aTL(ipl),eTL(ipl),ITL(ipl)/pi*180,stdList(il));
            Stde(il,ipl)=std(de);
            Me=Stde(il,ipl)^2*NCE(ipl);
            Ae(il,ipl)=sqrt(2/pi*Me);
        end
    end

end

%% Plot
figure;

set(gcf,'Position',[400,100,700/4/0.618*2,700],'color','w');
row=4;
col=2;

BottomRetainWidth=0.05;
LeftRetainWidth=0.09;
Height=0.23;
Width=0.40;

markersize=10;
linewidth=1;
fontsize=15;

subplot(row,col,1);
h1=semilogy(stdList,AI(:,1),'k.-','markersize',markersize,'linewidth',linewidth);hold all;
h2=semilogy(stdList,Ae(:,1),'k.--','markersize',markersize,'linewidth',linewidth);
xxlim=get(gca,'xlim');
semilogy(xxlim,[AI_data(1) AI_data(1)],'r-');
semilogy(xxlim,[Ae_data(1) Ae_data(1)],'r--');
set(gca,'xticklabel',[]);
set(gca,'position',[LeftRetainWidth BottomRetainWidth+3*Height Width Height]);
% yylim=get(gca,'ylim');
yylim=[1e-5 1e-2];
ylim(yylim);
set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));
ylabel('$|A_I|(|A_e|)$','fontsize',fontsize,'Interpreter','latex');
title([fnameList{1},'&2004UP10'],'fontsize',fontsize/5*4);
legend([h1 h2],{'$|A_e|$','$|A_I|$'},'fontsize',fontsize,'location','best','Interpreter','latex');

subplot(row,col,3);
h3=semilogy(stdList,StdI(:,1),'k.-','markersize',markersize,'linewidth',linewidth);hold all;
h4=semilogy(stdList,Stde(:,1),'k.--','markersize',markersize,'linewidth',linewidth);
xxlim=get(gca,'xlim');
plot(xxlim,[StdI_data(1) StdI_data(1)],'r-');
plot(xxlim,[Stde_data(1) Stde_data(1)],'r--');
yylim=[1e-6 1e-4];
ylim(yylim);
% set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
set(gca,'position',[LeftRetainWidth BottomRetainWidth+2*Height Width Height]);
% yylim=get(gca,'ylim');
set(gca,'yTickLabel',{'10^{-6}','10^{-5}',''});
xlabel('$\sigma_f~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
ylabel('$s_I(s_e)$','fontsize',fontsize,'Interpreter','latex');
legend([h3 h4],{'$s_e$','$s_I$'},'fontsize',fontsize,'location','best','Interpreter','latex');

% subplot(2,2,3);
% plot(stdList,0,'w-');hold all;
% xxlim=get(gca,'xlim');
% plot(xxlim,[NCE NCE],'k-');
% plot(xxlim,[NCE_data NCE_data],'r-');

subplot(row,col,2);
semilogy(stdList,AI(:,2),'k.-','markersize',markersize,'linewidth',linewidth);hold all;
semilogy(stdList,Ae(:,2),'k.--','markersize',markersize,'linewidth',linewidth);
xxlim=get(gca,'xlim');
semilogy(xxlim,[AI_data(2) AI_data(2)],'r-');
semilogy(xxlim,[Ae_data(2) Ae_data(2)],'r--');
set(gca,'xticklabel',[]);
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+3*Height Width Height]);
% yylim=get(gca,'ylim');
yylim=[1e-5 1e-2];
ylim(yylim);
title([fnameList{2},'&2004UP10'],'fontsize',fontsize/5*4);

set(gca,'yTick',power(10,log10(yylim(1))+1:log10(yylim(2))));
ylabel('$|A_I|(|A_e|)$','fontsize',fontsize,'Interpreter','latex');

subplot(row,col,4);
semilogy(stdList,StdI(:,2),'k.-','markersize',markersize,'linewidth',linewidth);hold all;
semilogy(stdList,Stde(:,2),'k.--','markersize',markersize,'linewidth',linewidth);
xxlim=get(gca,'xlim');
plot(xxlim,[StdI_data(2) StdI_data(2)],'r-');
plot(xxlim,[Stde_data(2) Stde_data(2)],'r--');
set(gca,'position',[2*LeftRetainWidth+Width BottomRetainWidth+2*Height Width Height]);
yylim=get(gca,'ylim');
set(gca,'yTick',power(10,log10(yylim(1)):log10(yylim(2))-1));
xlabel('$\sigma_f~\rm(DEG)$','fontsize',fontsize,'Interpreter','latex');
ylabel('$s_I(s_e)$','fontsize',fontsize,'Interpreter','latex');
