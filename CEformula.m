clear;
% aT=30;aP=40;eP=1-aT/aP;IP=3;
mP=6.56e-9;
aT=30;aP=40;eP=0.26;%IP=2;
R0=3.5*aP*(mP/3)^(1/3);
% R0=1/1.5e8;
% alpha=10;

IP=1:1:30;
PF=zeros(length(IP),1);
for i=1:length(IP)
    PF(i)=Fun_CollisionProb(aT,aP,eP,IP(i),R0);
end

figure;
semilogy(IP,PF*1e9,'k-');

