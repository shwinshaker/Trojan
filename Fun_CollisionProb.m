function [PF]=Fun_CollisionProb(aTL,aPL,ePL,IPL)

mu=(2*pi)^2*1^2;
mP=6.56e-9;

PF=zeros(length(aTL),1);
for i=1:length(aTL)
    aT=aTL(i);
    aP=aPL(i);
    eP=ePL(i);
    IP=IPL(i);
   
    R0=3.5*aP*(mP/3)^(1/3);
    
    A=aP/aT;
    cotAlpha2=(A^2*eP^2*(2*A-1)-(A-1)^2)/A^2/(1-eP^2);
    
    % vP=(mu*(2/aT-1/aP))^(1/2);
    % vT=(mu*(1/aT))^(1/2);
    % tP=2*pi/sqrt(mu/aP^3);
    % tT=2*pi/sqrt(mu/aT^3);
    
    U2=mu/aT^2*(aP*(2*A-1)/A^2+aT-2*(aP*aT*(1-eP^2))^(1/2)*cosd(IP));
    U=U2^(1/2);
    
    % TestFactor=R0/pi/aT/sind(IP)/cotAlpha2^(1/2)*(cotAlpha2+sind(IP)^2)^(1/2);
    % P1=aT/4/aT^2/(1-eP^2)^(1/2)*R0*U/((mu/aT)^(1/2)*A^(1/2)*(1-eP^2)^(1/2)*(cotAlpha2+sind(IP)^2)^(1/2));
    
    % Tisser=1/aP+2*sqrt(aP*(1-eP^2))*cosd(IP);
    % U=sqrt(3-Tisser);
    
    
    % TestFactor*pi*R0*sqrt(vP^2+vT^2-2*vP*vT*cosd(alpha))/vP/vT/sind(alpha)/tP/tT*1e18/(1.5e8)^2
    
    PF(i)=1/2/pi^2*R0^2*U*aT/(sind(IP)*cotAlpha2^(1/2)*aT^2*aP^2*(1-eP^2)^(1/2));

end
% PF*1e18
% 
% %%??????????
% Fun_NCE(aP,aT,eP,0,IP/180*pi,0,1)/1e9*1e18/(1.5e8)^2

end
