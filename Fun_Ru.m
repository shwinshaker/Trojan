function [RuL]=Fun_Ru(aPL,aTL,ePL,eTL,IPL,ITL)

%% Following is the intial method from Mixed_MassHeapDmin
% %% determine Ru
% %[x,y]=solve('x^2/aP^2+y^2/aP^2/(1-eP^2)=1','((x-aP*eP)*cos(IP-IT))^2+y^2=aT^2','x','y');
% [x,y]=solve('(x/cos(IP)+aP*eP)^2/aP^2+y^2/aP^2/(1-eP^2)=1','(x/cos(IT))^2+y^2=aT^2','x','y');
% yy=eval(y);
% for i=1:length(yy)
%     if isreal(yy(i)) && yy(i) > 0
%         WW=yy(i)*2;
%     end
% end
% 
% %% four types of averge
% % HH=((aP*(1-eP)*sin(IP)+aT*sin(IT))^2+(aP*(1-eP)*sin(IP)-aT*sin(IT))^2)^(1/2);
% % HH=(abs((aP*(1-eP)*sin(IP)+aT*sin(IT))*(aP*(1-eP)*sin(IP)-aT*sin(IT))))^(1/2);
% HH=(abs(aP*(1-eP)*sin(IP)+aT*sin(IT))+abs(aP*(1-eP)*sin(IP)-aT*sin(IT)))/2;
% % HH=2/(1/abs(aP*(1-eP)*sin(IP)+aT*sin(IT))+1/abs(aP*(1-eP)*sin(IP)-aT*sin(IT)));
% 
% 
% LL=abs(aP*(1-eP)*cos(IP)-aT*cos(IT));
% RR=(HH*LL*WW)^(1/3);
% TP=1/(1/aP^3)^(1/2);
% Ru=RR/(1e9/2/TP)^(1/2);
% disp('Ru:');

RuL=zeros(length(aPL),1);
for iL=1:length(aPL)
    aP=aPL(iL);
    aT=aTL(iL);
    eP=ePL(iL);
    eT=eTL(iL);
    IP=IPL(iL);
    IT=ITL(iL);
    %% determine Ru
    %% ??IP???IP?? eP???????
    %% ??????Trojan??????IT=20??????
    %% Trojan?????????????????????? ???????
%     [~,y]=solve('(x/cos(IP+IT)+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
%         'x^2+y^2=aT^2','x','y');
    %% ???????????

%     [~,y]=solve('(x/cos(IP)+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
%         '(x/cos(IT))^2+y^2=aT^2','x','y');
%     [~,y]=solve('(x+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
%         'x^2+y^2/(1-eT^2)=aT^2','x','y');
        [~,y]=solve('(x+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
        'x^2+y^2=aT^2','x','y');

    %% ???????????
%     [~,y]=solve('(x/cos(IP-IT)+aP*eP)^2/aP^2+y^2/aP^2/(1-eP^2)=1','x^2+y^2=aT^2','x','y');
    %xx=eval(x);
    WW=2*pi*aT*2; %% if no solution
    %WW=Inf;
    yy=eval(y);
    for i=1:length(yy)
        if isreal(yy(i)) && yy(i) > 0
            %WW=yy(i)*2*2;
            %xxx=xx(i);
            WW=asin(yy(i)/aT)*aT*2*2;
        end
    end
    
%     HH=2*(abs(aP*(1-eP)*sin(IP)+aT*sin(IT))+...,
%         abs(aP*(1-eP)*sin(IP)-aT*sin(IT)))/2;
%     HH=2*sqrt(abs(aP*(1-eP)*sin(IP)+aT*sin(IT))*...,
%         abs(aP*(1-eP)*sin(IP)-aT*sin(IT)));
    %% ???Trojan??
    HH=2*sqrt(abs(aP*(1-eP)*sin(IP+IT))*...,
        abs(aP*(1-eP)*sin(IP-IT)));

%     HH=2*sqrt(abs(aP*(1-eP)*IP+aT*IT)*...,
%         abs(aP*(1-eP)*IP-aT*IT));
%     HH=2*(abs(aP*(1-eP)*IP+aT*IT)+...,
%         abs(aP*(1-eP)*IP-aT*IT))/2;
%     HH=2*abs(aP*(1-eP)*sin(IP)-aT*sin(IT));
%     HH=2*(aP*(1-eP)+aT)/2*sin(abs(IP-IT));
%     HH=2*(abs(aP*(1-eP)+aT))/2*abs(IP-IT);

%     %% ????
%     HH1=2*sqrt(aP^2*(1-eP)^2+aT^2-2*aP*(1-eP)*aT*cos(IP-IT));
%     HH2=2*sqrt(aP^2*(1-eP)^2+aT^2-2*aP*(1-eP)*aT*cos(IP+IT));
%     HH=sqrt(HH1*HH2);

%     LL=abs(aP*(1-eP)*cos(IP)-aT*cos(IT));
%     LL=sqrt(abs(aP*(1-eP)-aT*(1-eT))*abs(aP*(1-eP)-aT*(1+eT)))*2;
    %% ?????????? ??CE??????????
    %% ??????????? ???????????CE?
    %% ???????????????????????? ????
    LL=abs(aP*(1-eP)-aT*(1-eT))+abs(aP*(1-eP)-aT*(1+eT));
%     LL=abs(aP*(1-eP)*cos(IP)-aT*(1-eT)*cos(IT))+abs(aP*(1-eP)*cos(IP)-aT*(1+eT)*cos(IT));

    %% LL?????????? eP?????????? ????????
    %% ???????????? ???eP?????? ??????????? ??????
    %% ??????? ????????? ?????????????????
    %% ????????????ep?????? ???????
%     disp(aP*(1-eP));
%     disp(aT*(1-eT));
%     disp(aT*(1+eT));
%     if aP*(1-eP)<aT*(1-eT)
%         LL=(aT-aP*(1-eP))*2;
%     elseif aP*(1-eP)>=aT*(1-eT) && aP*(1-eP)<aT*(1+eT)
%         LL=(aT*(1+eT)-aT*(1-eT))*2;
%     elseif aP*(1-eP)>=aT*(1+eT)
%         LL=inf;
%     end
    

    %% eP?????????????????? 
    %% ????????????????????
    %     disp(xxx);
%     LL=abs(max(aT,aP*(1-eP))-xxx)/2;
    
    RR=(HH*LL*WW/(4*pi/3))^(1/3);
    
    TP=1/(1/aP^3)^(1/2);
    Ru=RR/(1e9/2/TP)^(1/2);
    
    %     disp('Ru:');
    %     disp(Ru)
    RuL(iL)=Ru;
end

