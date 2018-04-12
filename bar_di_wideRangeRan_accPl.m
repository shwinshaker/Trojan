% CE times - total inc
tag='data';
%tag='repeat';
if ~strcmp(tag,'repeat')
clear;
Dirname={'pl'};

pl_name={'2001KN77';'1999CE119';'1998WV31';'2001FU172';'1993SB'; ...,
    '2001UO18';'1994TB';'1999CM158';'2000YH2';'2002CE251';'2000FB8';'1996TP66';'1998HQ151';...,
    '2001RX143';'1995QY9';'1996SZ4';...,
    '2000GE147';'2001KD77';'2001FR185';'1998HK151';...,
    '2002CW224';'2001KY76';'2000CK105';'1995HM5'};

npl_name={'2001KN77';'1999CE119';'1995QY9';'1998HK151';'2000FB8'};
npl20_name={'1999CE119','2001KN77'};
npl40_name={'1999CE119','2001KN77','1999CE119&2006RJ103','2001FU172&2006RJ103'};

ran_npl_name={'1999CE119_10k','1999CE119_100k','1999CE119_200k','1999CE119_300k','1999CE119_400k','1999CE119_500k','2001KN77_100k','1999CE119_1M',...,
    '1999CE119_100k_30000Fchk','2004UP10_100k_WideRange'};
Ndata=0;
for ii=1:length(Dirname)
    filename=eval([Dirname{ii},'_name']);
    Ndata=Ndata+length(filename);
end
% CE_times=cell(Ndata,1);
% Totdi=cell(Ndata,1);
% Maxdi=cell(Ndata,1);
% Meanabsdi=cell(Ndata,1);
% Varabsdi=cell(Ndata,1);
% Meandi=cell(Ndata,1);
% Vardi=cell(Ndata,1);
% MaxTotdi=cell(Ndata,1);
% name=cell(Ndata,1);
maxa=zeros(Ndata,1);
mina=zeros(Ndata,1);
maxe=zeros(Ndata,1);
mine=zeros(Ndata,1);
maxi=zeros(Ndata,1);
mini=zeros(Ndata,1);


Dir='ServerMount';
fontsize=15;

di_acc=[];
plel_acc=[];
tpel_acc=[];

for ii=1:length(Dirname)
    disp(Dirname{ii});
    filename=eval([Dirname{ii},'_name']);
    if strcmp(Dirname{ii},'pl')
        ffname='RealPlutinos';
        di_name='di_record_perturb';
        CE_name='CE_record';
    elseif strcmp(Dirname{ii},'ran_npl')
        ffname='RanPlutinos';
        di_name='di_fit_perturb';
        CE_name='ran_record';
    elseif intersect(Dirname{ii},['npl' 'npl20' 'npl40'])
        ffname='RealPlutinosNpl';
        di_name='di_record_perturb';
        CE_name='CE_record';
    end
    for i=1:length(filename)
        name{Ndata+i}=[filename{i},'_',Dirname{ii}];
        if strcmp(Dirname{ii},'ran_npl')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'npl20')
            fname=[filename{i},'_1Gyr_20pl'];
        elseif strcmp(Dirname{ii},'npl40')
            fname=[filename{i},'_1Gyr_40pl'];
        else
            fname=[filename{i},'_1Gyr'];
        end
        disp(filename{i});
        di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_name,'.txt'));
        di_acc=[di_acc;di_record_perturb];
        
        plel=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/plel.txt'));
        tpel=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/tpel.txt'));

        maxa(i)=max(plel(:,2));
        mina(i)=min(plel(:,2));
        maxe(i)=max(plel(:,3));
        mine(i)=min(plel(:,3));        
        maxi(i)=max(plel(:,4));
        mini(i)=min(plel(:,4));        
        plel_acc=[plel_acc;plel];
        tpel_acc=[tpel_acc;tpel];
    end
end
end

plel_ran=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/2004UP10_100k_WideRange2/AE_record_pl.txt'));
tpel_ran=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/2004UP10_100k_WideRange2/AE_record_tp.txt'));

plel_ran2=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/1999CE119_100k/AE_record_pl.txt'));
tpel_ran2=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/1999CE119_100k/AE_record_tp.txt'));

figure(1);
set(gcf,'Position',[400,100,700,400],'color','w');

subplot(2,2,1);
plot(39.5,1,'w');hold all;%axis equal;
ylim([0,0.5]);
xlim([39 40]);
plot(plel_ran(:,1),plel_ran(:,2),'y.');
for i=1:Ndata
    rectangle('position',[mina(i),mine(i),maxa(i)-mina(i),maxe(i)-mine(i)]);
end
plot(plel_acc(:,2),plel_acc(:,3),'k.');
plot(plel_ran2(:,1),plel_ran2(:,2),'r.');
hold off;

subplot(2,2,2);
plot(39.5,1,'w');hold all;%axis equal;
%ylim([0,40]);
xlim([39 40]);
plot(plel_ran(:,1),plel_ran(:,3),'y.');
for i=1:Ndata
    rectangle('position',[mina(i),mini(i),maxa(i)-mina(i),maxi(i)-mini(i)]);
end
plot(plel_acc(:,2),plel_acc(:,4),'k.');
plot(plel_ran2(:,1),plel_ran2(:,3),'r.');
hold off;

subplot(2,2,3);
plot(tpel_ran(:,1),tpel_ran(:,2),'y.');hold all;
plot(tpel_acc(:,2),tpel_acc(:,3),'k.');
plot(tpel_ran2(:,1),tpel_ran2(:,2),'r.');
hold off;

subplot(2,2,4);
plot(tpel_ran(:,1),tpel_ran(:,3),'y.');hold all;
plot(tpel_acc(:,2),tpel_acc(:,4),'k.');
plot(tpel_ran2(:,1),tpel_ran2(:,3),'r.');
hold off;


di_ran=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/2004UP10_100k_WideRange2/di_fit_perturb.txt'));

di_ran_100k=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/1999CE119_100k/di_fit_perturb.txt'));
%di_ran_100k=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RanPlutinos/2004UP10_100k_WideRange/di_fit_perturb.txt'));

di_npl=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/RealPlutinosNpl/1999CE119_1Gyr_40pl/di_record_perturb.txt'));

        
yylim=0.1;
N=50;

figure(2);
set(gcf,'Position',[400,100,700,500],'color','w');

Maxdi=max(abs(di_acc));
di_norm=di_acc/Maxdi;
di_norm_0=di_norm(di_norm~=0);
    
% H=kstest2(di_norm_0(di_norm_0>0),-di_norm_0(di_norm_0<0));
% disp(titlename{isub});
% disp(H);
    
Minlog=-log(min(abs(di_norm_0)));

dix=0:N;
dix=dix';
dlog=Minlog/N;
dix=exp(-dix*dlog);
dix=[-dix;dix];
dix=sort(dix);

countx=histcounts(di_norm,dix)'/length(di_acc);

di_norm_ran=di_ran/Maxdi;
countx_ran=histcounts(di_norm_ran,dix)'/length(di_ran);

di_norm_npl=di_npl/Maxdi;
countx_npl=histcounts(di_norm_npl,dix)'/length(di_npl);

di_norm_100k=di_ran_100k/Maxdi;
countx_100k=histcounts(di_norm_100k,dix)'/length(di_ran_100k);

di_norm_tail=di_record_perturb/Maxdi;
countx_tail=histcounts(di_norm_tail,dix)'/length(di_record_perturb);

dix=(dix(2:end)+dix(1:end-1))/2;

plotx_p=dix(dix>=0);
ploty_p=countx(dix>=0);
plotx_n=-dix(dix<0);
ploty_n=countx(dix<0);
    
ploty_p_ran=countx_ran(dix>=0);
ploty_n_ran=countx_ran(dix<0);

ploty_p_npl=countx_npl(dix>=0);
ploty_n_npl=countx_npl(dix<0);

ploty_p_100k=countx_100k(dix>=0);
ploty_n_100k=countx_100k(dix<0);

ploty_p_tail=countx_tail(dix>=0);
ploty_n_tail=countx_tail(dix<0);

semilogx(0,0,'w');hold all;

semilogx(plotx_p,ploty_p,'r.-','linewidth',1.3,'markersize',10);
semilogx(plotx_n,ploty_n,'b.-','linewidth',1.3,'markersize',10);
semilogx(plotx_p,ploty_p_ran,'r.--','linewidth',1.3,'markersize',10);
semilogx(plotx_n,ploty_n_ran,'b.--','linewidth',1.3,'markersize',10);
semilogx(plotx_p,ploty_p_npl,'r.-.','linewidth',1.3,'markersize',10);
semilogx(plotx_n,ploty_n_npl,'b.-.','linewidth',1.3,'markersize',10);
semilogx(plotx_p,ploty_p_100k,'r+--','linewidth',1.3,'markersize',10);
semilogx(plotx_n,ploty_n_100k,'b+--','linewidth',1.3,'markersize',10);
semilogx(plotx_p,ploty_p_tail,'ro--','linewidth',1.3,'markersize',10);
semilogx(plotx_n,ploty_n_tail,'bo--','linewidth',1.3,'markersize',10);

xlim([1e-8 10]);
ylim([0 yylim]);

%     absMeandi=abs(mean(di_norm));
%     Vardi=var(di_norm);
%     plot([absMeandi absMeandi],[0 yylim],'k--');
% 
%     Meanabsdi=mean(abs(di_norm));
%     Varabsdi=var(abs(di_norm));
%     plot([Meanabsdi Meanabsdi],[0 yylim],'k--','linewidth',2.0);
hold off;

set(text(2e-8,8/9*yylim,['$$N_{acc}= ',num2str(length(di_acc)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');
set(text(2e-8,7/9*yylim,['$$N_{ran}= ',num2str(length(di_ran)),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');

%title(titlename,'fontsize',fontsize);
set(gca,'xTick',power(10,-8:1:2));

xlabel('$\Delta i~\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Proportion$','fontsize',fontsize,'Interpreter','latex');

