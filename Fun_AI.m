function [AI,Ae,Aa]=Fun_AI(aPL,aTL,ePL,eTL,IPL,ITL,nuP)

NCE=Fun_NCE_PlAsTp(aPL,aTL,ePL,eTL,IPL/180*pi,ITL/180*pi,nuP);

%[dinc,decc]=Fun_diDstb_theo(100000,1,aP,eP,IP/pi*180,aT,eT,IT/pi*180);
N=100000;
[dinc,de,da] = Fun_diDstb_theo_PlAsTp(N,nuP,aPL,ePL,IPL,aTL,eTL,ITL);
dinc=dinc/180*pi;

MI=(std(dinc)*nuP.^(2/3)).^2.*NCE;
Me=(std(de)*nuP.^(2/3)).^2.*NCE;
Ma=(std(da)*nuP.^(2/3)).^2.*NCE;

AI=sqrt(2/pi*MI)/pi*180;
Ae=sqrt(2/pi*Me);
Aa=sqrt(2/pi*Ma);



