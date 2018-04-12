
clear;

f=@(x)2*abs(x).*(x.^2+1).^(-2);
F=@(x)x.^2./(1+x.^2);

N=1e6;
x=rand(N,1);
xx=rand(N,1);
y=(xx-1/2)./abs(xx-1/2).*(1./x-1).^(1/2);

% % probability density theo and exper
% y=abs(y);
% 
% Minlog=-log(min(y));
% Maxlog=log(max(y));
% 
% Nbin=100;
% yy=0:Nbin;
% yy=yy';
% yy=exp(-Minlog+yy*(Maxlog+Minlog)/Nbin);
% yy=sort(yy);
% 
% county=histcounts(y,yy)';
% ploty=(yy(1:end-1)+yy(2:end))/2;
% lenBin=yy(2:end)-yy(1:end-1);
% 
% semilogx(ploty,county/N./lenBin,'k.');hold all;
% semilogx(ploty,f(ploty),'r-');
% semilogx(ploty,F(ploty),'b-');

% simulated random walk
plot(cumsum(y),'k-');hold all;
plot(cummax(abs(y)),'r-');
plot(sqrt(1:N-1),'r-');


disp(max(abs(y)));

%% solve max value (find at least once beyond the value)
% syms a;
% vpasolve(1-F(a)==1/N,a);

% ymax=(N-1)^(1/2);