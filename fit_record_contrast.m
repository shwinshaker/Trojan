% % plot(CErelativeposition(:,1),CErelativeposition(:,2),'bo');
% % hold on;
% % plot(CErelativepositionfit(:,1),CErelativepositionfit(:,2),'r.');
% % axis([-2 2 -2 2]);axis square;
% % hold on;
% % plot(0,0,'k.')
% % theta=0:0.1:(2*pi+0.2);
% % R= (0.26)^(1/2)*3.5;
% % plot(R*cos(theta),R*sin(theta),'k');
% % a=4.749858746825433E-003;
% % e=363.367710358107;
% % b=a*(e^2-1)^(1/2);
% % plot(a*sec(theta),b*tan(theta),'r');
% 
% 
% figure;
% plot(p(1),p(2),'ko',pin(1),pin(2),'ko',pout(1),pout(2),'ko');
% axis([-0.3 0.3 -0.3 0.3]);axis square;
% hold on;plot(0,0,'+');
% text(p(1),p(2),'CE');
% 
%%3D trail and threshold sphere
figure;
plot3(tp(:,1)-pl(:,1),tp(:,2)-pl(:,2),tp(:,3)-pl(:,3),'ko');
axis([-0.15 0.15 -0.15 0.15 -0.15 0.15]);
axis square;grid on;grid minor;box on;
hold on;
plot3(0,0,0,'b+');hold on;
%plot3(pCE(1),pCE(2),pCE(3),'r+');
r=0.1798;[x y z]=sphere;mesh(r*x,r*y,r*z);alpha(0.2);

% %3D fit&record relative position
% figure;plot3(fit(:,1),fit(:,2),fit(:,3),'k-','markersize',0.1);
% hold on;plot3(record(:,1),record(:,2),record(:,3),'r-','markersize',0.1);
% box on;axis square;grid on;grid minor;
% hold on;plot3(0,0,0,'b+');hold on;
% r=0.1798;[x y z]=sphere;mesh(r*x,r*y,r*z);alpha(0.5);

%distance plot
% r=0.1798;
% figure;for i=1:length(fit1);plot(i/length(fit1),norm(fit1(i,:))/r,'k.','markersize',0.1);hold on;end;
% hold on;
% for j=1:length(record1);plot(j/length(record1),norm(record1(j,:))/r,'r.','markersize',0.1);hold on;end;
% axis([0 1 0 1]);axis square;grid on;grid minor;
% for i=1:length(fit4);plot(i/length(fit4),norm(fit4(i,:))/r,'k.','markersize',0.1);hold on;end;
% hold on;
% for j=1:length(record4);plot(j/length(record4),norm(record4(j,:))/r,'r.','markersize',0.1);hold on;end;


%3D fit&record relative position
%figure;
% plot3(tp(1,1),tp(1,2),tp(1,3),'ko','markersize',0.2);
% plot3(pl(1,1),pl(1,2),pl(1,3),'ro','markersize',0.2);
% axis([-30 0 -4 -3 0.2 0.3]);
% box on;axis square;grid on;grid minor;
% hold on;
% for i=1:length(tp)
% plot3(tp(i,1),tp(i,2),tp(i,3),'k.','markersize',0.1);
% axis([-40 -20 -10 -3 0.2 0.4]);
% box on;axis square;grid on;grid minor;
% hold on;
% plot3(pl(i,1),pl(i,2),pl(i,3),'r.','markersize',0.1);
% pause(0.1);
% end
% 

% hold on;plot3(0,0,0,'b+');hold on;
% r=0.1798;[x y z]=sphere;mesh(r*x,r*y,r*z);alpha(0.5);