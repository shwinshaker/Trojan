t=0:0.03:4*pi;
tp=0:0.02:4*pi;
tp1=zeros(1,length(tp));
y=sin(t);x=cos(t);
e=0.27;
yp=4/3*(1-e^2)^(1/2)*sin(tp);
xp=4/3*cos(tp)+4/3*e;
for ip1=1:length(tp)
    if ip1<=length(t)
        tp1(ip1)=t(ip1);
    else
        tp1(ip1)=t(ip1-length(t));
    end
end
xp1=xp.*cos(tp1)+yp.*sin(tp1);
yp1=-xp.*sin(tp1)+yp.*cos(tp1);
plot(x,y,'k');axis equal;axis([-2 2 -2 2]);hold all;
plot(xp,yp,'k-');
plot(xp1,yp1,'b-');
plot(x(1),y(1),'.r','markersize',50);
%hp=line('color','b','linestyle','-.','xdata',0,'ydata',1,'markersize',40,'erasemode','xor');
hp=plot(x(1),y(1),'.b','markersize',40,'erasemode','xor');
hp1=plot(xp(1),yp(1),'.k','markersize',40,'erasemode','xor');
ip=1;np=length(tp);
while 1
    set(hp,'xdata',xp1(ip),'ydata',yp1(ip));
    set(hp1,'xdata',xp(ip),'ydata',yp(ip));
    drawnow;
    pause(0.01);
    ip=ip+1;
    if ip>np
        ip=1;
    end
end
hold off;