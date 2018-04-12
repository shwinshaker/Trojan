function [ASim] = Fun_totEff_sim(nuP,aP,eP,IP,aT,eT,IT)

Ru=Fun_Ru(aP,aT,eP,eT,IP,IT);
% gm=3.5;
mP=6.56e-9;

%% Calculate ChiN and ChiM
ChiN=(3.5*aP/Ru)^2*(mP/3)^(2/3);
A=(3-aT/aP)/2;
B=(2-aT/aP)^(1/2);
ChiM=(aT/Ru)^2*(2/(A-B*cos(max(IP,IT))))*mP^2;
MSim=ChiM*nuP.^2.*(-1+log(ChiN)+log(nuP)/3*2);
ASim=sqrt(2/pi*MSim);

end