% xhin=-29.195;xhout=-29.086;
% yhin=-3.849;yhout=-5.988;
% vxhin=0.0000772;vxhout=0.000279;
% vyhin=-0.00354575;vyhout=-0.00351166;
% md=0.1798;
% V=5.25e-4;
% alpha=142.6;
% beta=52.8;
% a=0:1:90;
% b=0;

xhin=17.633;xhout=18.462;
yhin=-25.049;yhout=-24.173;
vxhin=0.00240;vxhout=0.00233;
vyhin=-0.00245;vyhout=-0.00254;
md=0.1798;
V=7.5469e-04;
alpha=118.1749;
beta=77.3579;
a=0:1:90;
b=0;

xhtin=xhin+md*cosd(alpha+a);
yhtin=yhin+md*sind(alpha+a);
vxhtin=vxhin+V*cosd(beta-b);
vyhtin=vyhin+V*sind(beta-b);

xhtout=xhout+md*cosd(alpha-a);
yhtout=yhout+md*sind(alpha-a);
vxhtout=vxhout+V*cosd(beta+b);
vyhtout=vyhout+V*sind(beta+b);

miu=(2*pi)^2/365.25^2;

rin=(xhtin.^2+yhtin.^2).^(1/2);
vin=(vxhtin.^2+vyhtin.^2).^(1/2);
hin=abs(xhtin.*vyhtin-yhtin.*vxhtin);
ain=1./(2./rin-vin.^2/miu);
ein=sqrt(1-hin.^2./(miu.*ain));

rout=(xhtout.^2+yhtout.^2).^(1/2);
vout=(vxhtout.^2+vyhtout.^2).^(1/2);
hout=abs(xhtout.*vyhtout-yhtout.*vxhtout);
aout=1./(2./rout-vout.^2./miu);
eout=sqrt(1-hout.^2./(miu.*aout));

de=eout-ein;

figure;
plot(a,de,'k-');
set(gca,'XDir','reverse');%?X????