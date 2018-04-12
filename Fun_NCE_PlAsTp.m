function [NCE]=Fun_NCE_PlAsTp(aPL,aTL,ePL,eTL,IPL,ITL,nuP)

gm=3.5;
mP=6.56e-9;

Ru=Fun_Ru(aPL,aTL,ePL,eTL,IPL,ITL);

NCE=(gm*aTL./Ru).^2*(mP/3)^(2/3)*nuP.^(2/3);