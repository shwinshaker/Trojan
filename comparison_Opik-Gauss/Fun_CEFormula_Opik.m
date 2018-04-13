function [di,de,da,xb,yb,zb,sinPhi,cosPhi,x0,y0,z0,Ux,Uy,Uz] = Fun_CEFormula_Opik(aP,OP,mP,at,et,it,Ot,wt,ft)
% Given pre-encounter elements, calculate the post-encounter elements
% Planet: circular planar
% Particle: no limits
% based on Valsecchi 2015 & 2003 & Carusi 1990

% time from ft
if ft==0
    error('f should not be 0, i.e. particle should not be at perihelion');
end
if at*(1-et)>aP || at*(1+et)<aP
    error('q too large or Q too small');
end
t=0;

disp(' ');
disp('-------------------------');
disp('      --- Opik ---       ')
disp('-------------------------');
fprintf('Pre-at: %.4f AU\n',at);
fprintf('Pre-et: %.4f\n',et);
fprintf('Pre-it: %.4f DEG\n',it);
disp('-------------------------');
fprintf('Pre-OP: %.4f DEG\n',OP);
fprintf('Pre-Ot: %.4f DEG \n',Ot);
fprintf('Pre-wt: %.4f DEG\n',wt);
fprintf('Pre-ft: %.4f DEG\n',ft);
disp('-------------------------');
dth=3.5*(mP/3).^(1/3);
fprintf('Rth/aP: %.4f\n', dth);

%% 
cosit=cosd(it);sinit=sind(it);
cosft=cosd(ft);sinft=sind(ft); % cosft=sqrt(1-sinft^2)^(1/2);
coswt=cosd(wt);sinwt=sind(wt);

coswft=coswt*cosft-sinwt*sinft;
sinwft=sinwt*cosft+coswt*sinft;
r=at.*(1-et.^2)./(1+et.*cosft);
beta=asind(sinwft.*sinit);
% lambda=Ot+2*atand(sinwft.*cosit./(coswft+sqrt(1-sinwft.^2.*sinit.^2)));
% deltaLambda=lambda-OP;
%% Valsecchi 2015
% angy=sinwft.*cosit;
% angx=coswft+sqrt(1-sinwft.^2.*sinit.^2);
% angnorm=sqrt(angx.^2+angy.^2);
% cosang=angx./angnorm;
% sinang=angy./angnorm;
% cos2ang=cosang.^2-sinang.^2;
% sin2ang=2*cosang.*sinang;
% cosOtMnsOP=cosd(Ot-OP);
% sinOtMnsOP=sind(Ot-OP);
% cosdeltaLambda=cosOtMnsOP.*cos2ang-sinOtMnsOP.*sin2ang;
% sindeltaLambda=sinOtMnsOP.*cos2ang+cosOtMnsOP.*sin2ang;
%% My formula based on Napier's rules
%%% 在非常小的半径的奇点上表现更好
angy=sinwft.*cosit;
angx=coswft;
angnorm=sqrt(angx.^2+angy.^2);
cosang=angx./angnorm;
sinang=angy./angnorm;
cosOtMnsOP=cosd(Ot-OP);
sinOtMnsOP=sind(Ot-OP);
cosdeltaLambda=cosOtMnsOP.*cosang-sinOtMnsOP.*sinang;
sindeltaLambda=sinOtMnsOP.*cosang+cosOtMnsOP.*sinang;

%% to XYZ
% x=r./aP.*cosd(deltaLambda).*cosd(beta)-1;
% y=r./aP.*sind(deltaLambda).*cosd(beta);
x0=r./aP.*cosdeltaLambda.*cosd(beta)-1;
y0=r./aP.*sindeltaLambda.*cosd(beta);
z0=r./aP.*sind(beta);
fprintf('Enter position: (%.6f,%.6f,%.6f)\n', x0,y0,z0);
% disp(x);

%% to opik variable
U=sqrt(3-aP./at-2*sqrt(at.*(1-et.^2)./aP).*cosit);
cosTheta=(1-U.^2-aP./at)/2./U;
sinTheta=sqrt(2-aP./at-at.*(1-et.^2).*cosit.^2./aP)./U;
sign=(sinft>0)*2-1;
sign2=(coswft>0)*2-1;
sinPhi=sign.*sqrt(2-aP./at-at.*(1-et.^2)./aP)./U./sinTheta;
cosPhi=sign2.*sqrt(at.*(1-et.^2)./aP).*sinit./U./sinTheta;
Ux=U.*sinTheta.*sinPhi;
Uy=U.*cosTheta;
Uz=U.*sinTheta.*cosPhi;
fprintf('relative velocity: %.4f\n', U);
% fprintf('f sign: %i\n',sign);
% fprintf('cos phi: %.4f\n',cosPhi);
% fprintf('sin phi: %.4f\n',sinPhi);
fprintf('Pre-Phi: %.4f\n',angle(cosPhi+sinPhi*1i)/pi*180);
% fprintf('cos theta: %.4f\n',cosTheta);
% fprintf('sin theta: %.4f\n',sinTheta);
fprintf('Pre-Theta: %.4f\n',angle(cosTheta+sinTheta*1i)/pi*180);

disp('-------------------------')
%% to b-plane
itm=(x0.*sinPhi+z0.*cosPhi).*sinTheta+y0.*cosTheta;
tb=t-itm./U;
fprintf('tb-t*: %.4f\n',tb-t); 
%%% tb 正负应该无所谓吧 也不能确定随便选的初始点是CE前还是后啊 只要CE了就算嘛

xb=x0-itm.*sinTheta.*sinPhi;
yb=y0-itm.*cosTheta;
zb=z0-itm.*sinTheta.*cosPhi;
fprintf('Pre-CE position: (%.6f,%.6f,%.6f)\n', xb,yb,zb);
% minimum distance on asymptote
db=sqrt(xb.^2+yb.^2+zb.^2);
fprintf('Pre-minimum distance (asymptotic): %.4f\n', db);
if db>dth
    fprintf('R0/Rth: %.4f\n', db/dth);
    error('R0>Rth: too far!')
end

xi=xb.*cosPhi-zb.*sinPhi;
zeta=(xb.*sinPhi+zb.*cosPhi).*cosTheta-yb.*sinTheta;

disp('-------------------------')
%% swing-by
c=mP./U.^2;
c2=c.^2;
b2=xi.^2+zeta.^2;
b2Pc2=b2+c2;
b2Mc2=b2-c2;
b2Mc2SinTheta=b2Mc2.*sinTheta;
cZeta=c.*zeta;
cXi=c.*xi;
TwoCZetaCosTheta=2*cZeta.*cosTheta;
TwoCZetaSinTheta=2*cZeta.*sinTheta;
TwoCXiCosTheta=2*cXi.*cosTheta;
TwoCXiSinTheta=2*cXi.*sinTheta;
Us=U;
cosThetas=(b2Mc2.*cosTheta+TwoCZetaSinTheta)./b2Pc2;
sinThetas=sqrt((b2Mc2SinTheta-TwoCZetaCosTheta).^2+4*c2.*xi.^2)./b2Pc2;
cosPhis=((b2Mc2SinTheta-TwoCZetaCosTheta).*cosPhi+TwoCXiSinTheta)./b2Pc2./sinThetas;
sinPhis=((b2Mc2SinTheta-TwoCZetaCosTheta).*sinPhi-TwoCXiCosTheta)./b2Pc2./sinThetas;
xis=xi.*sinTheta./sinThetas;
zetas=(b2Mc2SinTheta.*zeta-2*b2.*c.*cosTheta)./b2Pc2./sinThetas;
tbs=tb;
Uxs=Us.*sinThetas.*sinPhis;
Uys=Us.*cosThetas;
Uzs=Us.*sinThetas.*cosPhis;
fprintf('Post-Thetas: %.4f\n',angle(cosThetas+sinThetas*1i)/pi*180);
% fprintf('sinThetas: %.4f\n',sinThetas);
fprintf('Post-Phi: %.4f\n',angle(cosPhis+sinPhis*1i)/pi*180);
% fprintf('cosPhis: %.4f\n',cosPhis);
% fprintf('sinPhis: %.4f\n',sinPhis);
if cosThetas>1
    fprintf('>>>>>>>>>>>>>>>>> Larger than 1 Force to be 1!\n');
    fprintf('cosThetas: %.8f\n', cosThetas);
    cosThetas=1;
%     error('larger than 1');
end
if sinThetas>1
    fprintf('sinThetas: %.8f\n', sinThetas);
    error('larger than 1');
end
if cosPhis>1
    fprintf('>>>>>>>>>>>>>>>>> Larger than 1 Force to be 1!\n');
    fprintf('cosPhis: %.8f\n', cosPhis);
    cosPhis=1;
%     error('larger than 1');
end
if sinPhis>1
    fprintf('sinPhis: %.8f\n', sinPhis);
    error('larger than 1');
end

%% reverse
% fluctuation 在 a e i 上出现的还不一样: a比较少 e次之 i最多 
% 猜想是角度计算上的问题，因为三者用到的角度不一样
ats=aP./(1-Us.^2-2*Us.*cosThetas);
fprintf('Post-at: %.4f\n',ats);
ets=Us.*sqrt((Us+2*cosThetas).^2+sinThetas.^2.*sinPhis.^2.*(1-Us.^2-2*Us.*cosThetas));
fprintf('Post-et: %.4f\n',ets);
cosits=(1+Us.*cosThetas)./sqrt(1+2*Us.*cosThetas+Us.^2.*(1-sinThetas.^2.*sinPhis.^2));
Uzs=Us.*sinThetas.*cosPhis;
Uys=Us.*cosThetas;
sinits=sqrt(Uzs.^2./(Uzs.^2+(1+Uys).^2)); % Carusi 1990
%% Choose which its to use: Use "angle" function for accuracy
its_flag='cos';
switch its_flag
    case 'cos' 
        its=acosd(cosits);
    case 'sin'
        its=asind(sinits);
    case 'cos+sin'
        its=angle(cosits+sinits*1i)/pi*180;
end
fprintf('Post-it: %.4f\n',its);
%% 

xbs=zetas.*cosThetas.*sinPhis+xis.*cosPhis;
zbs=zetas.*cosThetas.*cosPhis-xis.*sinPhis;
ybs=-zetas.*sinThetas;
fprintf('Post-CE position: (%.6f,%.6f,%.6f)\n', xbs,ybs,zbs);
dbs=sqrt(xbs.^2+ybs.^2+zbs.^2);
fprintf('Post-minimum distance (asymptotic): %.4f\n', dbs);

xmin=(xb+xbs)/2;
ymin=(yb+ybs)/2;
zmin=(zb+zbs)/2;
fprintf('P0: (%.6f,%.6f,%.6f)\n', xmin,ymin,zmin);
dmin=sqrt(xmin.^2+ymin.^2+zmin.^2);
fprintf('Average-minimum distance (asymptotic): %.4f\n', dmin);
% These two should be perpendicular
fprintf('x*v: %.8f\n',sum([Ux,Uy,Uz].*[xb,yb,zb]));
fprintf('xs*vs: %.8f\n',sum([Uxs,Uys,Uzs].*[xbs,ybs,zbs]));
% This should not be perpendicular
fprintf('xs*vs: %.8f\n',sum([Uxs,Uys,Uzs].*[xmin,ymin,zmin]));
% Compare three distances
fprintf('d0/db %.6f\n',dmin/db);
fprintf('d0/ds %.6f\n',dmin/dbs);
fprintf('r norm %.6f\n',norm([xbs;ybs;zbs]-[xb;yb;zb]));

disp('-------------------------');
di=its-it;
fprintf('di %.6f\n',di);
de=ets-et;
fprintf('de %.6f\n',de);
da=ats-at;
fprintf('da %.6f\n',da);





