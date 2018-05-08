% function [dinc,de,thetap,phip,thetav,phiv,phivpvt,lambda,theta,phi,f] = Fun_diDstb_theo(N,nuP,aP,eP,IP,aT,eT,IT,std0)
function [dinc,de,da] = Fun_diDstb_theo(nuP,aP,eP,IP,aT,eT,IT)

%% azm
[theta,phi,~,~,~,~,phivpvt,lambda,f,ki,N]=...,
    Fun_azmDstb_theo(aP,eP,IP,aT,IT);

coswf=cosd(lambda);

cosphivpvt=cosd(phivpvt);

cosf=cosd(f);
sinf=sind(f);
cosE=(eT+cosf)./(1+eT*cosf);

%% Coeff
%std0=std(plel(:,4));
mPluto=6.56e-9;
mP=nuP*mPluto;
gm=3.5;
Rth=gm*aP*(mP/3)^(1/3);

A=(3-aT/aP)/2;
B=(2-aT/aP)^(1/2);
% dinc0=mP*aT/Rth*(2/(1-eT^2))^(1/2);
dinc0=mP*aT/Rth*(2*(1-eT^2))^(1/2);
de0=mP*aT/Rth*(2*(1-eT^2))^(1/2);
da0=aT*mP*aT/Rth*2*(2/(1-eT^2))^(1/2);

%% linear Dstb of R
syms xx;
F=@(x)(x^2); %% cdf
G=matlabFunction(finverse(F(xx)));

a=0;b=1;
y=a+rand(N,1)*(b-a);
gmR=G(y);

%%% resulting di Dstb
%%% mainR=(1./gmR.^2-1).^(1/2);
mainR=1./gmR;
% dinc=dinc0*coswf.*sind(theta)./(A-B*cosphivpvt).^(1/2).*(1./gmR.^2-1).^(1/2)/pi*180;
rhoi=coswf./(1+eT.*cosf);
dinc=dinc0*rhoi.*sind(theta)./(A-B*cosphivpvt).^(1/2).*mainR/pi*180;
zetae=sinf.*cosd(phi)+cosf.*sind(phi)+cosE.*sind(phi);
de=de0*zetae.*cosd(theta)./(A-B*cosphivpvt).^(1/2).*mainR;
zetaa=eT*(sinf.*cosd(phi)+cosf.*sind(phi))+sind(phi);
da=da0*zetaa.*cosd(theta)./(A-B*cosphivpvt).^(1/2).*mainR;

end

