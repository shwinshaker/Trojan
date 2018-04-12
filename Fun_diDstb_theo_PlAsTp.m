% function [dinc,de,thetap,phip,thetav,phiv,phivpvt,lambda,theta,phi,f] = Fun_diDstb_theo(N,nuP,aP,eP,IP,aT,eT,IT,std0)
function [dinc,de,da] = Fun_diDstb_theo_PlAsTp(N,nuP,aP,eP,IP,aT,eT,IT)

%% azm
std0=abs(max(IT,IP))^(1/3)+1;

[theta,phi,~,~,~,~,phivpvt,lambda,f]=...,
    Fun_azmDstb_theo(N,aP,eP,IP,aT,IT,std0);
% [theta,phi,~,~,~,~,phivpvt,lambda,f]=...,
%     Fun_azmDstb_theo(N,aT,eP,IT,aP,IP,std0);

coswf=cosd(lambda);

cosphivpvt=cosd(phivpvt);

cosf=cosd(f);
sinf=sind(f);
cosE=(eP+cosf)./(1+eP*cosf);

%% Coeff
%std0=std(plel(:,4));
mPluto=6.56e-9;
mP=nuP*mPluto;
gm=3.5;
Rth=gm*aT*(mP/3)^(1/3);

A=(3-aT/aP)/2;
B=(2-aT/aP)^(1/2);
dinc0=mP*aT*(aT/aP)^(1/2)/Rth*(2/(1-eP^2))^(1/2);
de0=mP*(aP*aT)^(1/2)/Rth*(2*(1-eP^2))^(1/2);
da0=mP*aP^2*(aT/aP)^(1/2)/Rth*2*(2/(1-eT^2))^(1/2);

%% linear Dstb of R
syms xx;
F=@(x)(x^2); %% cdf
G=matlabFunction(finverse(F(xx)));

a=0;b=1;
y=a+rand(N,1)*(b-a);
gmR=G(y);

%%% resulting di Dstb
dinc=dinc0*coswf.*sind(theta)./(A-B*cosphivpvt).^(1/2).*(1./gmR.^2-1).^(1/2)/pi*180;
zetae=sinf.*cosd(phi)+cosf.*sind(phi)+cosE.*sind(phi);
de=de0*zetae.*cosd(theta)./(A-B*cosphivpvt).^(1/2).*(1./gmR.^2-1).^(1/2);
zetaa=eP*(sinf.*cosd(phi)+cosf.*sind(phi))+sind(phi);
da=da0*zetaa.*cosd(theta)./(A-B*cosphivpvt).^(1/2).*(1./gmR.^2-1).^(1/2);

end

