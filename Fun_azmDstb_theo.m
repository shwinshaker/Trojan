function [theta,phi,thetap,phip,thetav,phiv,phivpvt,lambda,f,ki,N] = ...,
    Fun_azmDstb_theo(aP,eP,IP,aT,IT)

N=1000000;

%% free parameter
% std0=2.2*abs(IP)^0.5;
std0=abs(IP)^(0.6)+3.5;
% std0=3.5;

% std1=std0;
% std0=abs(IP)^(1.0)+2;
% std0=2.5*abs(IP)^0.45;
% std0=abs(IP)^(0.7)+2;

fprintf('std: %.4f\n',std0);

%% Theo
% meanp=IP; %% ?????1???
% %stdp=8.3568/2.4;
% stdp=std0;%plincstd;
% % f(thetap)=exp(-(plotx-meanvp).^2/2/stdvp^2)/(2*pi*stdvp^2)^(1/2);
% 

% stdpos=std0;%7.2041/5;
% stdneg=stdpos;

%% random number generation

% thetap=[normrnd(IP+IT,std0,[round(N/4) 1]);...,
%     normrnd(-IP-IT,std0,[round(N/4) 1]);...,
%     normrnd(IP-IT,std0,[round(N/4) 1]);...,
%     normrnd(-IP+IT,std0,[N-3*round(N/4) 1])];

% phippos=normrnd(meanpos,stdpos,[round(N/2) 1]);
% phipneg=normrnd(meanneg,stdneg,[N-round(N/2) 1]);

kappa=1/circ_ang2rad(std0)^2;
thetap=[circ_vmrnd(circ_ang2rad(IP+IT),kappa,round(N/4));...,
        circ_vmrnd(circ_ang2rad(-IP-IT),kappa,round(N/4));...,
        circ_vmrnd(circ_ang2rad(IP-IT),kappa,round(N/4));...,
        circ_vmrnd(circ_ang2rad(-IP+IT),kappa,N-3*round(N/4))];
thetap=circ_rad2ang(thetap);

%% phip
% deltaang_max=Fun_deltaang(aP,eP,IP,aT);
% deltaang_min=Fun_deltaang(aP,eP,0,aT);
% meanpos_max=90-deltaang_max;
% meanneg_max=90+deltaang_max;
% meanpos_min=90-deltaang_min;
% meanneg_min=90+deltaang_min;

% kappa=1/circ_ang2rad(std1)^2;
% phip=[circ_vmrnd(circ_ang2rad(meanpos_max),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanpos_min),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanneg_max),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanneg_min),kappa,N-3*round(N/4))];
% phip=circ_rad2ang(phip);

% deltaang=[Fun_deltaang(aP,eP,IP+IT,aT);...,
%           Fun_deltaang(aP,eP,IP-IT,aT)];
% deltaang=deltaang/cosd(IP/2);

% meanphi=[90-deltaang(1);90+deltaang(1);90-deltaang(2);90+deltaang(2)];
% phip=[circ_vmrnd(circ_ang2rad(meanphi(1)),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanphi(2)),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanphi(3)),kappa,round(N/4));...,
%         circ_vmrnd(circ_ang2rad(meanphi(4)),kappa,N-3*round(N/4))];
% phip=circ_rad2ang(phip);

deltaang=Fun_deltaang(aP,eP,IP,aT);
meanpos=90-deltaang;
meanneg=90+deltaang;
phip=[circ_vmrnd(circ_ang2rad(meanpos),kappa,round(N/2));...,
      circ_vmrnd(circ_ang2rad(meanneg),kappa,N-round(N/2))];
phip=circ_rad2ang(phip);


%% calc
if exist('i','var')
    clear i;
end
        
Vp=(2/aT-1/aP)^(1/2);
Vt=(1/aT)^(1/2);

SH=Vp*cosd(thetap);
PH=Vp*sind(thetap);
HT=(SH.^2+Vt^2-2*SH*Vt.*sind(phip)).^(1/2);
% thetav=atand(PH./HT);

vx=Vp*cosd(thetap).*cosd(phip);
vy=Vp*cosd(thetap).*sind(phip)-Vt;
vz=Vp*sind(thetap);
% thetav=atand(vz./(vx.^2+vy.^2).^(1/2));
thetav=circ_rad2ang(angle((vx.^2+vy.^2).^(1/2)+vz*1i));
% phiv=atand(vy./vx);
% phiv(vx<0 & vy<0)=phiv(vx<0 & vy<0)+180;
% phiv(vx<0 & vy>0)=phiv(vx<0 & vy>0)+180;
% % phiv(vx>0 & vy<0)=phiv(vx>0 & vy<0)+360;
phiv=circ_rad2ang(angle(vx+vy*1i));


%% Dstb of phiv
PT=(PH.^2+HT.^2).^(1/2);
cosphivpvt=(Vp^2+Vt^2-PT.^2)/(2*Vp*Vt);
phivpvt=acosd(cosphivpvt);
% phivpvt=acosd(cosd(thetap).*sind(phip));

%% Dstb of theta
% f(ki)=Const;
ki=rand(N,1)*360;
theta=asind(sind(ki).*cosd(thetav));

%% Dstb of phi
zeta=abs(atand(sind(thetav)./cotd(ki)));
phivmod=mod(phiv,360);
% figure;
% histogram(zeta);
phi=zeros(length(phiv),1);
% phi(phivmod<180)=phiv(phivmod<180)-(90+zeta(phivmod<180));
% phi(phivmod>180)=phiv(phivmod>180)+(90+zeta(phivmod>180));
phi(ki>270)=phivmod(ki>270)-(90+zeta(ki>270));
phi(ki>180 & ki<270)=phivmod(ki>180 & ki<270)+(90+zeta(ki>180 & ki<270));
phi(ki>90 & ki<180)=phivmod(ki>90 & ki<180)+(90-zeta(ki>90 & ki<180));
phi(ki<90)=phivmod(ki<90)-(90-zeta(ki<90));
phi=mod(phi,360);

% phi=mod(phivmod-(90-zeta),360);

%% Dstb of coswf
std=8*abs(IT-IP-4)+25;
kappa=1/circ_ang2rad(std)^2;
lambda_pos=circ_vmrnd(circ_ang2rad(0),kappa,round(N/2));
lambda_pos=circ_rad2ang(lambda_pos);
lambda_neg=circ_vmrnd(circ_ang2rad(180),kappa,N-round(N/2));
lambda_neg=circ_rad2ang(lambda_neg);

lambda=[lambda_pos;lambda_neg];

% lambda=rand(N,1)*360-180;
% lambda=rand(N,1)*360;

%% Dstb of f
f=rand(N,1)*360-180;

