clear;
syms x;
F=@(x)(x^2); %% cdf
G=matlabFunction(finverse(F(x)));

N=100000;
a=0;b=1;
y=a+rand(N,1)*(b-a);
x=G(y);

Nbin=50;
xbin=(0:Nbin)'/Nbin;
plotx=xbin(2:end)-1/N/2;
counts=histcounts(x,xbin);
ploty=counts/N;
plot(plotx,ploty,'k.-');hold on;
plot(plotx,0.04*plotx,'r-');