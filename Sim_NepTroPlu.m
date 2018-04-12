aNep=30.134;aTro=30.099;aPlu=39.583;
ePlu=0.274;eTro=0.025;
nNep=(1/aNep^3)^(1/2);
nTro=(1/aTro^3)^(1/2);
nPlu=(1/aPlu^3)^(1/2);
tNep=0:nNep:4*pi;
tTro=0:nTro:4*pi;
tTro=tTro+60/180*pi;
tPlu=0:nPlu:4*pi;
tPlu=tPlu-15/180*pi;


yNep=aNep*sin(tNep);xNep=aNep*cos(tNep);
yTro=aTro*(1-eTro^2)^(1/2)*sin(tTro);xTro=aTro*(cos(tTro)+eTro);
yPlu=aPlu*(1-ePlu^2)^(1/2)*sin(tPlu);xPlu=aPlu*(cos(tPlu)+ePlu);

tNep1=zeros(1,length(tPlu));
for i=1:length(tPlu)
    if i<=length(tNep)
        tNep1(i)=tNep(i);
    else
        tNep1(i)=tNep(i-length(tNep));
    end
end
xPlu1=xPlu.*cos(tNep1)+yPlu.*sin(tNep1);
yPlu1=-xPlu.*sin(tNep1)+yPlu.*cos(tNep1);

tNep2=zeros(1,length(tTro));
for i=1:length(tTro)
    if i<=length(tNep)
        tNep2(i)=tNep(i);
    else
        tNep2(i)=tNep(i-length(tNep));
    end
end
xTro1=xTro.*cos(tNep2)+yTro.*sin(tNep2);
yTro1=-xTro.*sin(tNep2)+yTro.*cos(tNep2);

plot(0,0,'w+');hold all;axis equal;axis square;axis([-60 60 -60 60]);
plot(xNep,yNep,'k-');
plot(xTro,yTro,'r-');
plot(xPlu,yPlu,'b-');
plot(xPlu1,yPlu1,'b--');
plot(xTro1,yTro1,'r--');

hNep=plot(xNep(1),yNep(1),'.k','markersize',50);
hTro=plot(xTro(1),yTro(1),'.r','markersize',40);
hTro1=plot(xTro1(1),yTro1(1),'.m','markersize',40);
hPlu=plot(xPlu(1),yPlu(1),'.b','markersize',40);
hPlu1=plot(xPlu1(1),yPlu1(1),'.c','markersize',40);

ip=1;ip1=1;
it=1;it1=1;
iN=1;
while 1
    set(hPlu,'xdata',xPlu(ip),'ydata',yPlu(ip));
    set(hPlu1,'xdata',xPlu1(ip1),'ydata',yPlu1(ip1));
    set(hTro,'xdata',xTro(it),'ydata',yTro(it));
    set(hTro1,'xdata',xTro1(it1),'ydata',yTro1(it1));
    set(hNep,'xdata',xNep(iN),'ydata',yNep(iN));
    drawnow;
    pause(0.001);
    ip=ip+1;
    ip1=ip1+1;
    it=it+1;
    it1=it1+1;
    iN=iN+1;
    if ip>length(xPlu)
        ip=1;
    end
    if ip1>length(xPlu1)
        ip1=1;
    end
    if it>length(xTro)
        it=1;
    end
    if it1>length(xTro1)
        it1=1;
    end
    if iN>length(xNep)
        iN=1;
    end
end
% hold off;