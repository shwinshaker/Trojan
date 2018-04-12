%% Test N-D distribution
N=1000000;
%% linear Dstb of R
syms xx;
F=@(x)(x^2); %% cdf
G=matlabFunction(finverse(F(xx)));

a=0;b=1;
y=a+rand(N,1)*(b-a);
gmR=G(y);

%% uniform Dstb of theta
a=-pi/2;b=pi/2;
theta=a+rand(N,1)*(b-a);

%% uniform Dstb of lambda
a=0;b=2*pi;
lambda=a+rand(N,1)*(b-a);

z=1e-6*sin(theta).*cos(lambda).*(1./gmR.^2-1).^(1/2);

figure;

Nbin=400;
Minlog=-log(min(abs(z)));

dix=(0:Nbin)';
dlog=Minlog/Nbin;
dix=exp(-dix*dlog);
dix=[-dix;dix];
dix=sort(dix);

plotx=(dix(2:end)+dix(1:end-1))/2;
dixlen=dix(2:end)-dix(1:end-1);

countx=histcounts(z,dix)'/N;

semilogx(plotx(plotx>0),countx(plotx>0),'k.-','linewidth',1.3,'markersize',10);hold all;
semilogx(-plotx(plotx<0),countx(plotx<0),'k.--','linewidth',1.3,'markersize',10);

f=@(z)1/2./(1+z.^2).^(3/2);
semilogx(plotx(plotx>0),f(plotx(plotx>0)/1e-6).*dixlen(plotx>0)/1e-6,'r.-','linewidth',1.3,'markersize',10);
semilogx(-plotx(plotx<0),f(plotx(plotx<0)/1e-6).*dixlen(plotx<0)/1e-6,'b.-','linewidth',1.3,'markersize',10);


