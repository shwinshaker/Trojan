function [OPL,OtL,wtL,ftL] = Fun_CEFormula_GenRand(aP,mP,at,et,it,Nran)
%% Generate random CE simulation

% time from ft
if at*(1-et)>aP || at*(1+et)<aP
    error('q too large or Q too small');
end

dth=3.5.*(mP/3).^(1/3);  %% aP is cancalled by R0 scale
cosit=cosd(it);sinit=sind(it);
atmul1mnsesqr=at*(1-et^2);
aPdivat=aP/at;

disp('-------------------------');
disp('     --- Ran Gen ---     ')
disp('-------------------------');
fprintf('aP: %.4f AU\n',aP);
fprintf('at: %.4f AU\n',at);
fprintf('et: %.4f \n',et);
fprintf('it: %.4f DEG\n',it);

% ran=rand(Nran,3)*360;
% OPL=ran(:,1);OtL=ran(:,2);wtL=ran(:,3);
OPL=zeros(Nran,1);OtL=zeros(Nran,1);
wtL=zeros(Nran,1);ftL=zeros(Nran,1);
searchRange=20; %DEG
% ftL=OPL-mod(OtL+wtL,360)+(rand(Nran,1)*2-1)*searchRange;
% ixL=zeros(Nran,1);
ixpp=0;

% for ix=1:Nran
while ixpp < Nran
%     OP=OPL(ix);Ot=OtL(ix);wt=wtL(ix);ft=ftL(ix);
    %% Ran Gen
    OP=rand()*360;Ot=rand()*360;wt=rand()*360;
    ft=OP-mod(Ot+wt,360)+(rand()*2-1)*searchRange;
    
    %% 
    cosft=cosd(ft);sinft=sind(ft); % cosft=sqrt(1-sinft^2)^(1/2);
    coswt=cosd(wt);sinwt=sind(wt);
    coswft=coswt*cosft-sinwt*sinft;
    sinwft=sinwt*cosft+coswt*sinft;
    r=atmul1mnsesqr/(1+et*cosft);
    rdivaP=r/aP;
    sinbeta=sinwft*sinit;
    cosbeta=sqrt(1-sinbeta^2);
    
%     %% Valsecchi 2015
%     lambda=Ot+2*atand(sinwft*cosit/(coswft+sqrt(1-sinwft^2.*sinit^2)));
%     deltaLambda=lambda-OP;
%     cosdeltaLambda=cosd(deltaLambda);
%     sindeltaLambda=sind(deltaLambda);
    
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
    x=rdivaP*cosdeltaLambda*cosbeta-1;
    y=rdivaP*sindeltaLambda*cosbeta;
    z=rdivaP*sinbeta;
    
    %% to opik variable
    U=sqrt(3-aPdivat-2*sqrt(atmul1mnsesqr/aP)*cosit);
    cosTheta=(1-U^2-aPdivat)/2/U;
    sinTheta=sqrt(2-aPdivat-atmul1mnsesqr*cosit^2/aP)/U;
    sign=(sinft>0)*2-1;
    sign2=(coswft>0)*2-1;
    sinPhi=sign*sqrt(2-aPdivat-atmul1mnsesqr/aP)/U/sinTheta;
    cosPhi=sign2*sqrt(atmul1mnsesqr/aP)*sinit/U/sinTheta;
    
    %% to b-plane
    itm=(x*sinPhi+z*cosPhi)*sinTheta+y*cosTheta;
    xb=x-itm*sinTheta*sinPhi;
    yb=y-itm*cosTheta;
    zb=z-itm*sinTheta*cosPhi;
    db=sqrt(xb.^2+yb.^2+zb.^2);
    if db < dth
        ixpp=ixpp+1;
        % fprintf('Detected: %i\n', ix);
        fprintf('\nDetected #: %i\n', ixpp);
        fprintf('Minimum distance: %.4f\n', db/dth);
        OPL(ixpp)=OP;OtL(ixpp)=Ot;wtL(ixpp)=wt;ftL(ixpp)=ft;
        % ixL(ixpp)=ix;
    end
end

% ixL(ixL==0)=[];
% OPLL=OPL(ixL);
% OtLL=OtL(ixL);
% wtLL=wtL(ixL);
% ftLL=ftL(ixL);






