% %%tp1000plot
% 
% %%tp1000 L4L5 
% %%input: tpTrojanel & Nepel
% NL4=0;NL5=0;
% figure;
% for i=1:length(tpTrojanel)
%     angle=mod(tpTrojanel(i,5)+tpTrojanel(i,6)+tpTrojanel(i,7)-(Nepel(1000,5)+Nepel(1000,6)+Nepel(1000,7)),360);
%     if (angle<0) 
%         angle=angle+360;
%     end
%     if(tpTrojanel(i,2)<31 && tpTrojanel(i,2)>29)
%     if (angle>0 && angle<180)
%         NL4=NL4+1;
%         subplot(2,2,1);
%         plot(angle,tpTrojanel(i,4),'ok','markersize',tpTrojanel(i,3)*100);axis([0 180 0 30]);title('L4');
%         legend(strcat('N=',num2str(NL4)));
%         hold on;
% %         plot(angle,tpTrojanelinit(find(tpTrojanelinit(:,1)==tpTrojanel(i,1)),4)+sigmai(find(sigmai(:,1)==tpTrojanel(i,1)),2),'ob','markersize',tpTrojanel(i,3)*100);
%         subplot(2,2,3);
%         plot(tpTrojanel(i,3),tpTrojanel(i,4),'ok','markersize',(180-angle)/360*20);axis([0 0.4 0 30]);title('L4');
%         legend(strcat('N=',num2str(NL4)));
% %         plot(tpTrojanel(i,3),tpTrojanelinit(find(tpTrojanelinit(:,1)==tpTrojanel(i,1)),4)+sigmai(find(sigmai(:,1)==tpTrojanel(i,1)),2),'ob','markersize',angle/360*20);axis([0 0.4 0 30]);title('L4');
%         hold on;
%     else
%         NL5=NL5+1;
%         subplot(2,2,2);
%         plot(angle,tpTrojanel(i,4),'ok','markersize',tpTrojanel(i,3)*100);axis([180 360  0 30]);title('L5');
%         legend(strcat('N=',num2str(NL5)));
% %         plot(angle,tpTrojanelinit(find(tpTrojanelinit(:,1)==tpTrojanel(i,1)),4)+sigmai(find(sigmai(:,1)==tpTrojanel(i,1)),2),'ob','markersize',tpTrojanel(i,3)*100);
%         hold on;
%         subplot(2,2,4);
%         plot(tpTrojanel(i,3),tpTrojanel(i,4),'ok','markersize',(angle-180)/360*20);axis([0 0.4 0 30]);title('L5');
%         legend(strcat('N=',num2str(NL5)));
% %         plot(tpTrojanel(i,3),tpTrojanelinit(find(tpTrojanelinit(:,1)==tpTrojanel(i,1)),4)+sigmai(find(sigmai(:,1)==tpTrojanel(i,1)),2),'ob','markersize',(angle-180)/360*20);
%         hold on;
%     end
%     end
% end
% figure;
% subplot(3,1,1);plot(Nepel(:,1),Nepel(:,2));axis([0 max(Nepel(:,1)) 29 31]);
% subplot(3,1,2);plot(Nepel(:,1),Nepel(:,3));axis([0 max(Nepel(:,1)) 0 0.2]);
% subplot(3,1,3);plot(Nepel(:,1),Nepel(:,4));axis([0 max(Nepel(:,1)) 0 15]);
figure;
subplot(3,1,1);plot(el(:,1),el(:,2));axis([0 max(el(:,1)) 29 31]);
subplot(3,1,2);plot(el(:,1),el(:,3));axis([0 max(el(:,1)) 0 0.2]);
subplot(3,1,3);plot(el(:,1),el(:,4));axis([0 max(el(:,1)) 0 15]);

% %%sigmai
% for i=1:length(tpTrojanel)
% sigmai_clone(i,2)=sigmai(find(sigmai(:,1)==tpTrojanel(i,1)),2);
% sigmai_clone(i,1)=tpTrojanelinit(find(tpTrojanelinit(:,1)==tpTrojanel(i,1)),4);
% end
% sub1=sigmai_clone;sub2=sigmai_clone;
% sub1(sigmai_clone(:,2)<0,2)=0;
% sub2(sigmai_clone(:,2)>0,2)=0;
% subplot(2,1,1);
% h1=bar(sub1,'stack');ylim([-5 20]);set(h1(2),'facecolor','y');hold on;
% h2=bar(sub2,'stack');set(h2(2),'facecolor','m');
% subplot(2,1,2);
% bar(tpTrojanel(:,4),'r');ylim([-5 20]);
%Lamda-Lamda Libration Angle
%  figure;
%  plot(elOnlyNep(:,1),mod(elOnlyNep(:,5)+elOnlyNep(:,6)+elOnlyNep(:,7)-(NepOnlyel(:,5)+NepOnlyel(:,6)+NepOnlyel(:,7)),360),'k');
%  axis([0 max(elOnlyNep(:,1)) 0 360]);hold on;
%  plot(get(gca,'xlim'),[180 180],get(gca,'xlim'),[60 60],'r',get(gca,'xlim'),[300 300],'r');
%  
%  figure;
%  subplot(3,1,1);
%  plot(elOnlyNep(:,1),elOnlyNep(:,4));axis([0 max(elOnlyNep(:,1)) 0 20]);
%  subplot(3,1,2);
%  plot(elOnlyNep(:,1),elOnlyNep(:,2));axis([0 max(elOnlyNep(:,1)) 25 35]);
%  subplot(3,1,3);
%  plot(elOnlyNep(:,1),elOnlyNep(:,3));axis([0 max(elOnlyNep(:,1)) 0 0.1]);