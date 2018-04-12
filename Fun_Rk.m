function [Rk]=Fun_Rk(aPL,aTL,nuP)

% gm=3.5;
% VSN=(1/aT+mN/aT)^(1/2);
% VS=(1/aT)^(1/2);
mP=6.56e-9;
mN=5.17e-5;

Vcomp=(mN./aTL).^(1/2);
% Rk=16/15*mP*nuP./(((2./aTL-1./aPL).^(1/2)-(1./aTL).^(1/2)).*Vcomp);
Rk=2^(1/2)*mP*nuP./(((2./aTL-1./aPL).^(1/2)-(1./aTL).^(1/2)).*Vcomp);

%%% The previous notes for Rk
% % FR=@(nu)((-2*aT^2*(nu*mP).^2+...,
% %     2*(3^(1/3)*aP^2*aT^2*gm^2*(nu*mP).^(8/3)+3*aT^4*(nu*mP).^4).^(1/2)/3^(1/2)).^(1/2));
% %Vcomp=(1/aT)^(1/2);
% Vcomp=(mN/aT)^(1/2);
% % Vcomp=VSN-VS;
% FR=@(mu)16/15*mu/(((2/aT-1/aP)^(1/2)-(1/aT)^(1/2))*Vcomp);
% % FR=@(mu)mu./(((3^(1/3)/gm/aP*mu.^(2/3)).^(1/2)+(mN/aT)^(1/2)).^2);
% % FR=@(mu)mu./((VSN-VS+(VSN^2+3^(1/3)/gm/aP*mu.^(2/3)).^(1/2)).^2-VSN^2);
% FRS=@(mu)mu/mN*aT;
% % FRP=@(mP)(aP*exp(-3.9)/(4.5)^(1/2)*(mP/3)^(1/3)); %% 4.5Gyr not 1Gyr
% % lnx=PT1(2);
% % FRP=@(mP)(aP*gm*exp(-lnx/2)*(mP/3)^(1/3));
