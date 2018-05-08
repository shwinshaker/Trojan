function [di,de,da,Rr] = Fun_CEFormula_Gauss(aP,OP,mP,at,et,it,Ot,wt,xb,yb,zb,sinPhi,cosPhi,CEth)
% Given pre-encounter elements, calculate the element change
% Planet: circular planar
% Particle: no limits
% Rectilinear - Gauss
% Use [xb;yb;zb] from opik's b-plane as closest position

disp(' ');
disp('-------------------------');
disp('     --- Gaussian ---    ');
disp('-------------------------');

%% Calculate cosf sinf coswf sinwf at CE position according to Valsecchi 2015
cosf0=(at.*(1-et.^2)-aP.*(1+xb))./aP./et./(1+xb); %%% caution 1/et
cosf0=min(1,cosf0);  %% force to be less than 1
sign=(sinPhi>0)*2-1;
sinf0=sign.*sqrt(1-cosf0.^2);
cosE0=(et+cosf0)./(1+et.*cosf0);
sinwf0=aP.*(1+et.*cosf0).*zb./at./(1-et.^2)./sind(it);
sign=(cosPhi>0)*2-1;
coswf0=sign.*sqrt(1-sinwf0.^2);

%% rotate to the coordinate system co-moving with Trojan
%% So as to apply Gauss perturbation equation
P0=[cosd(OP) -sind(OP) 0;..., 
    sind(OP) cosd(OP) 0;...,
    0 0 1];
P1=[cosd(Ot) sind(Ot) 0;...,
    -sind(Ot) cosd(Ot) 0;...,
    0 0 1];
P2=[1 0 0;...,
    0 cosd(it) sind(it);...,
    0 -sind(it) cosd(it)];
P3=[cosd(wt) sind(wt) 0;...,
    -sind(wt) cosd(wt) 0;...,
    0 0 1];
P4=[coswf0 sinwf0 0;...,
    -sinwf0 coswf0 0;...,
    0 0 1];
% P0=[cosd(OP) sind(OP) 0;..., 
%     -sind(OP) cosd(OP) 0;...,
%     0 0 1];
% P1=[cosd(Ot) -sind(Ot) 0;...,
%     sind(Ot) cosd(Ot) 0;...,
%     0 0 1];
% P2=[1 0 0;...,
%     0 cosd(it) -sind(it);...,
%     0 sind(it) cosd(it)];
% P3=[cosd(wt) -sind(wt) 0;...,
%     sind(wt) cosd(wt) 0;...,
%     0 0 1];
% P4=[coswf0 -sinwf0 0;...,
%     sinwf0 coswf0 0;...,
%     0 0 1];
Pr=P4*P2*P1*P0*[xb;yb;zb];
x0=Pr(1);y0=Pr(2);z0=Pr(3);
% x0=xb;y0=yb;z0=zb;

%% Calculate the longitude and latitud of CE position
%%% On the plane of particle, whose element change is to be calculated
R0=sqrt(x0.^2+y0.^2+z0.^2)*aP;
fprintf('R0: %.4f AU\n',R0);
Rth=CEth*aP.*(mP/3).^(1/3);
fprintf('Rth: %.4f AU\n',Rth);
% if Rth<R0
%     error('R0>Rth: too far!')
% end
sinTheta0=z0./sqrt(x0.^2+y0.^2+z0.^2);
cosTheta0=sqrt(x0.^2+y0.^2)./sqrt(x0.^2+y0.^2+z0.^2);
sinPhi0=y0./sqrt(x0.^2+y0.^2);
cosPhi0=x0./sqrt(x0.^2+y0.^2);

% clear i;
fprintf('theta0: %.4f\n',angle(cosTheta0+sinTheta0*1i)/pi*180);
fprintf('phi0: %.4f\n', angle(cosPhi0+sinPhi0*1i)/pi*180);
fprintf('f0: %.4f\n', angle(cosf0+sinf0*1i)/pi*180);

%% calculate the angle between velocity vectors: alpha_v
mu=1; % mu 在这应该不影响，最后反正也是要消掉的
nt=sqrt(mu/at.^3);
% vxt=-nt.*at.*sind(ft)./sqrt(1-et.^2);
% vyt=nt.*at.*(et+cosd(ft))./sqrt(1-et.^2);
vxt=-nt.*at.*sinf0./sqrt(1-et.^2);
vyt=nt.*at.*(et+cosf0)./sqrt(1-et.^2);
vzt=0;
Vtr=[vxt;vyt;vzt];
nP=sqrt(mu/aP.^3);
vxP=-nP.*aP.*sind(OP);
vyP=nP.*aP.*cosd(OP);
vzP=0;
% P0=[cosd(OP) -sind(OP) 0;..., 
%     sind(OP) cosd(OP) 0;...,
%     0 0 1];
% P1=[cosd(Ot) sind(Ot) 0;...,
%     -sind(Ot) cosd(Ot) 0;...,
%     0 0 1];
% P2=[1 0 0;...,
%     0 cosd(it) sind(it);...,
%     0 -sind(it) cosd(it)];
% P3=[cosd(wt) sind(wt) 0;...,
%     -sind(wt) cosd(wt) 0;...,
%     0 0 1];
%%% 这里不要转P0了啊 因为速度已经是OP坐标系下的速度了啊少年
VPr=P3*P2*P1*[vxP;vyP;vzP];
cosAlpha=sum(VPr.*Vtr)./norm(VPr)./norm(Vtr);
fprintf('alpha_v: %.4f\n',acosd(cosAlpha));
% cosAlpha=cosd(acosd(cosAlpha)-30);

disp('-------------------------');
%% other arguments
% 这里要把p和t倒过来，因为p才是circular的，而t是有偏心率的
% A=(3-at/aP)/2;B=sqrt(2-at/aP);
A=(3-aP/at)/2;B=sqrt(2-aP/at);
% 因为p和t倒了，所以要除掉一个p乘一个t。因为默认提出去的是t，而实际应该是p
rltV=(A-B.*cosAlpha).^(1/2)/sqrt(aP)*sqrt(at);
% rltV=U;  % 这里直接用opik的U非常精确，但可能只是巧合
fprintf('rltV: %.4f\n', rltV);

%% headings
rt=at.*(1-et.^2)./(1+et.*cosf0);
% di0=-mP*at./Rth*(2./(1-et.^2)).^(1/2);
di0=-mP*rt./Rth*(2./(1-et.^2)).^(1/2);
de0=-mP*at./Rth*(2*(1-et.^2)).^(1/2);
da0=-at*mP*at./Rth*2*(2./(1-et.^2)).^(1/2);
% mainR=(1./(R0./Rth).^2-1).^(1/2);
mainR=1./(R0./Rth);

disp('-------------------------');
%% calc
di=di0.*coswf0.*sinTheta0./rltV.*mainR/pi*180;
fprintf('di1 %.6f\n',di);
zetae=sinf0.*cosPhi0+cosf0.*sinPhi0+cosE0.*sinPhi0;
de=de0*zetae.*cosTheta0./rltV.*mainR;
fprintf('de1 %.6f\n',de);
zetaa=et*(sinf0.*cosPhi0+cosf0.*sinPhi0)+sinPhi0;
da=da0*zetaa.*cosTheta0./rltV.*mainR;
fprintf('da1 %.6f\n',da);

if ~isreal(di) || ~isreal(de) || ~isreal(da)
    error('Complex!');
end

disp('-------------------------');
% fprintf('di/di1 %.4f\n',di/di1 )
% fprintf('dic/di1 %.4f\n',dic/di1)
% fprintf('dic2/di1 %.4f\n',dic2/di1)
% fprintf('de/de1 %.4f\n',de/de1 )
% fprintf('da/da1 %.4f\n',da/da1 )

Rr=R0/Rth;