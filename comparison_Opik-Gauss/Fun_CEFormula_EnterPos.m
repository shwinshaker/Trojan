function [x0,y0,z0] = Fun_CEFormula_EnterPos(aP,OP,at,et,it,Ot,wt,ft)
% Given pre-encounter elements, calculate the elements at the closest
% position
% Planet: circular planar
% Particle: no limits
% based on Valsecchi 2015

% time from ft
if at*(1-et)>aP || at*(1+et)<aP
    error('q too large or Q too small');
end

cosit=cosd(it);sinit=sind(it);
atmul1mnsesqr=at*(1-et^2);

%%
cosft=cosd(ft);sinft=sind(ft); % cosft=sqrt(1-sinft^2)^(1/2);
coswt=cosd(wt);sinwt=sind(wt);
coswft=coswt*cosft-sinwt*sinft;
sinwft=sinwt*cosft+coswt*sinft;
r=atmul1mnsesqr/(1+et*cosft);
rdivaP=r/aP;
sinbeta=sinwft*sinit;
cosbeta=sqrt(1-sinbeta^2);

%% Napier
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
x0=rdivaP*cosdeltaLambda*cosbeta-1;
y0=rdivaP*sindeltaLambda*cosbeta;
z0=rdivaP*sinbeta;