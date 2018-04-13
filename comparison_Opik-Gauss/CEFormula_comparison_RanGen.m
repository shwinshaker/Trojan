clear;
aP=30;
at=40;et=0.25;it=30;
mP=6e-9;% mP=3e-6;
Rth=3.5*aP.*(mP/3).^(1/3);

%% Rmin
mPluto=6e-9; % /Sun
RPluto=1188.3; %km
Rmin=((mP/mPluto)^(1/3)*RPluto)/(Rth*1.5e8);

%% generate angles for CE
%%% Load previous
if exist('CEFormula_cmp_RanGen_OPOtwtft.mat','file')
    load('CEFormula_cmp_RanGen_OPOtwtft.mat');
end
%% Get new
Nran=1000;
[OPL,OtL,wtL,ftL] = Fun_CEFormula_GenRand(aP,mP,at,et,it,Nran);
datmat2=[OPL,OtL,wtL,ftL];
%%% Merge
if exist('datmat','var')
    datmat=[datmat;datmat2];
else
    datmat=datmat2;
end
clear datmat2;
save('CEFormula_cmp_RanGen_OPOtwtft.mat','datmat');
%%
OPL=datmat(:,1);OtL=datmat(:,2);wtL=datmat(:,3);ftL=datmat(:,4);
len=length(OPL);

%% 0: Opik; 1: Gauss
di=zeros(len,1);di1=zeros(len,1);
de=zeros(len,1);de1=zeros(len,1);
da=zeros(len,1);da1=zeros(len,1);
Rr=zeros(len,1);

%% calculate element change
for ix=1:len
    OP=OPL(ix);Ot=OtL(ix);wt=wtL(ix);ft=ftL(ix);
    [di(ix),de(ix),da(ix),xb,yb,zb,sinPhi,cosPhi] = ...,
        Fun_CEFormula_Opik(aP,OP,mP,at,et,it,Ot,wt,ft);
    [di1(ix),de1(ix),da1(ix),Rr(ix)] = ...,
        Fun_CEFormula_Gauss(aP,OP,mP,at,et,it,Ot,wt,xb,yb,zb,sinPhi,cosPhi);
end

%% plot
figure;
set(gcf,'Position',[400,100,700,500],'color','w');
Nsub=3;
fontsize=20;

% before abs, check sign
if ~all(di.*di1) || ~all(de.*de1) || ~all(da.*da1)
    error('Sign not consistent!');
end
% abs for log
di=abs(di);di1=abs(di1);
de=abs(de);de1=abs(de1);
da=abs(da);da1=abs(da1);

%% Rr as argument
% 插值来把f写在上坐标轴把
subplot(Nsub,1,1)
semilogx(Rr,(di1-di)./di,'k.');hold all;
ylim([-2 2]);
xxlim=get(gca,'xlim');
plot(xxlim,[0 0],'k--');
% loglog(Rr,di1,'r.');
yylim=get(gca,'ylim');
plot([Rmin Rmin],yylim,'r-');
hold off;
xlabel('$R_0/R_{th}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\Delta i\mathrm{(DEG)}$','fontsize',fontsize,'Interpreter','latex');

subplot(Nsub,1,2)
semilogx(Rr,(de1-de)./de,'k.');hold all;
ylim([-2 2]);
xxlim=get(gca,'xlim');
plot(xxlim,[0 0],'k--');
% loglog(Rr,de1,'r.');
yylim=get(gca,'ylim');
plot([Rmin Rmin],yylim,'r-');
hold off;
xlabel('$R_0/R_{th}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\Delta e$','fontsize',fontsize,'Interpreter','latex');

subplot(Nsub,1,3)
plot(Rr,(da1-da)./da,'k.');hold all;
ylim([-2 2]);
xxlim=get(gca,'xlim');
plot(xxlim,[0 0],'k--');
% loglog(Rr,da1,'r.');
yylim=get(gca,'ylim');
plot([Rmin Rmin],yylim,'r-');
hold off;
xlabel('$R_0/R_{th}$','fontsize',fontsize,'Interpreter','latex');
ylabel('$\Delta a\mathrm{(AU)}$','fontsize',fontsize,'Interpreter','latex');

