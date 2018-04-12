clear;
aNep=30.134;aTro=30.099;aPlu=39.583;
ePlu=0.274;eTro=0.025;
nNep=(1/aNep^3)^(1/2);
nTro=(1/aTro^3)^(1/2);
nPlu=(1/aPlu^3)^(1/2);
tNep=0:nNep:4*pi;
tTro=0:nTro:4*pi;
tTro=tTro+60/180*pi;
tPlu=0:nPlu:4*pi;
tPlu=tPlu-20/180*pi;
% 
NtTro=length(tTro);
for i=NtTro+1:length(tPlu)
        tTro=[tTro tTro(i-NtTro)];
end

yNep=aNep*sin(tNep);xNep=aNep*cos(tNep);
yTro=aTro*(1-eTro^2)^(1/2)*sin(tTro);xTro=aTro*(cos(tTro)+eTro);
yPlu=aPlu*(1-ePlu^2)^(1/2)*sin(tPlu);xPlu=aPlu*(cos(tPlu)+ePlu);
yPlu=yPlu-yTro;
xPlu=xPlu-xTro;
ang=atan2(yTro,xTro);

ang1=zeros(1,length(tPlu));
for i=1:length(tPlu)
    if i<=length(ang)
        ang1(i)=ang(i);
    else
        ang1(i)=ang(i-length(ang));
    end
end
xPlu1=xPlu.*cos(ang1)+yPlu.*sin(ang1);
yPlu1=-xPlu.*sin(ang1)+yPlu.*cos(ang1);

xSun=-xTro;ySun=-yTro;
xSun1=xSun.*cos(ang1)+ySun.*sin(ang1);
ySun1=-xSun.*sin(ang1)+ySun.*cos(ang1);
% tNep2=zeros(1,length(tTro));
% for i=1:length(tTro)
%     if i<=length(tNep)
%         tNep2(i)=tNep(i);
%     else
%         tNep2(i)=tNep(i-length(tNep));
%     end
% end
% xTro1=xTro.*cos(tNep2)+yTro.*sin(tNep2);
% yTro1=-xTro.*sin(tNep2)+yTro.*cos(tNep2);
% xTro1=xTro.*cos(ang)+yTro.*sin(ang);
% yTro1=-xTro.*sin(ang)+yTro.*cos(ang);



%%%%%
plot(0,0,'w+');hold all;axis equal;axis square;axis([-80 80 -80 80]);
quiver(0,0,-xSun1(1),-ySun1(1),'r-','autoscale','on');
plot([0,xSun1(1)],[0,ySun1(1)],'r--');
plot(0,0,'k.','markersize',30);%% Trojan
% plot(xNep,yNep,'k-');
% plot(xTro,yTro,'k:');
plot(xSun1,ySun1,'k-','linewidth',1);
plot(xSun1(1),ySun1(1),'k.','markersize',50);%% Sun
% plot(xPlu,yPlu,'b-');
plot(xPlu1,yPlu1,'k-','linewidth',1); 
% plot(xTro1,yTro1,'k-','linewidth',1);
% plot(xTro1(1),yTro1(1),'k.','markersize',30);
plot(xPlu1(1),yPlu1(1),'k.','markersize',30); %% Plutino

fontsize=20;
xlabel('$x\rm (AU)$','fontsize',fontsize,'Interpreter','latex')
ylabel('$y\rm (AU)$','fontsize',fontsize,'Interpreter','latex')

% hNep=plot(xNep(1),yNep(1),'.k','markersize',50);
% hTro=plot(xTro(1),yTro(1),'.r','markersize',40);
% hTro1=plot(xTro(1),yTro(1),'.m','markersize',40);
% hPlu=plot(xPlu(1),yPlu(1),'.b','markersize',40);
% hPlu1=plot(xPlu1(1),yPlu1(1),'.c','markersize',40);
% 
% plot([0 xTro(1)],[0 yTro(1)],'r-');
% ip=1;ip1=1;
% it=1;it1=1;
% iN=1;
% while 1
%     set(hPlu,'xdata',xPlu(ip),'ydata',yPlu(ip));
%     set(hPlu1,'xdata',xPlu1(ip1),'ydata',yPlu1(ip1));
%     set(hTro,'xdata',xTro(it),'ydata',yTro(it));
% %     set(hTro1,'xdata',xTro1(it1),'ydata',yTro1(it1));
%     set(hNep,'xdata',xNep(iN),'ydata',yNep(iN));
%     drawnow;
%     pause(0.001);
%     ip=ip+1;
%     ip1=ip1+1;
%     it=it+1;
%     it1=it1+1;
%     iN=iN+1;
%     if ip>length(xPlu)
%         ip=1;
%     end
%     if ip1>length(xPlu1)
%         ip1=1;
%     end
%     if it>length(xTro)
%         it=1;
%     end
% %     if it1>length(xTro1)
% %         it1=1;
% %     end
%     if iN>length(xNep)
%         iN=1;
%     end
% end
% % hold off;