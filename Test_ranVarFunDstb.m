%% Test random variable function distribution
N=10000000;
x=rand(N,1)*2*pi;
y=cos(x);

xx=rand(N,1)*pi-pi/2;
yy=sin(xx);

z=y.*yy;

Nbin=1000;
xbin=(1:Nbin)'/Nbin*2-1;
plotx=(xbin(1:end-1)+xbin(2:end))/2;
binlen=xbin(2:end)-xbin(1:end-1);
countx=histcounts(z,xbin)'/N;

plot(plotx,countx./binlen);hold on;
f=@(x)1/pi./(1-x.^2).^(1/2);
plot(plotx,f(plotx),'r-');
