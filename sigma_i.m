
miuS=4*pi^2/365.25^2;
miu=miuS*6.5607561e-9*1000;
%x=3 y=4 z=5 vx=6 vy=7 vz=8
%tp2=trojan tp1=plutino
h=zeros(length(tp2ce(:,1)),3);
R=zeros(length(tp2ce(:,1)),3);V=zeros(length(tp2ce(:,1)),3);
RV=zeros(length(tp2ce(:,1)),3);RR=zeros(length(tp2ce(:,1)),3);
di=zeros(length(tp2ce(:,1)),1);
eh_record=zeros(length(tp2ce(:,1)),1);
%% const set !!!!!!!!!!!!!!
mode=2;%1=hyperbola 2=straight line;
md=1;%au %better to choose 3.5Rhill
dl=0.01;%au recommended length for each sub step / 0.001 is better
Nsubstep=100;
RRt=zeros(Nsubstep,3);RVt=zeros(Nsubstep,3);
dit=zeros(Nsubstep,1);
% rec_dcapm=0.1; % rad
for i=1:length(tp2ce) %length(tp2ce)
 R(i,:)=[tp2ce(i,3),tp2ce(i,4),tp2ce(i,5)];
 V(i,:)=[tp2ce(i,6),tp2ce(i,7),tp2ce(i,8)];
 RV(i,:)=[tp1ce(i,5)-tp2ce(i,6),tp1ce(i,6)-tp2ce(i,7),tp1ce(i,7)-tp2ce(i,8)];
 RR(i,:)=[tp1ce(i,2)-tp2ce(i,3),tp1ce(i,3)-tp2ce(i,4),tp1ce(i,4)-tp2ce(i,5)];
 a=tp2ceel(i,2);e=tp2ceel(i,3);inc=tp2ceel(i,4);capom=tp2ceel(i,5);beta=sqrt(1-e^2);
 sinwf=R(i,3)/(norm(R(i,:))*sind(inc));
 coswf=secd(capom)*(R(i,1)/norm(R(i,:))+sind(capom)*sinwf*cosd(inc));
 cosf=1/e*(a*(1-e^2)/norm(R(i,:))-1); 
 na=sqrt(miuS/a);
 h(i,:)=cross(R(i,:),V(i,:));
 if(mode==2)
% straight line approximation 
 det=sqrt(md^2-norm(RR(i,:))^2)/norm(RV(i,:));
%  Nsubstep=ceil(sqrt(md^2-norm(RR(i,:))^2)/dl);
 dt=det/Nsubstep;
 elseif(mode==1)
 %%hyperbola approximation
 %% Ellipse and Parabola Warning Preset
 [ialpha,ah,eh,inch,capomh,omegah,capmh]=xv2el(RR(i,1),RR(i,2),RR(i,3),RV(i,1),RV(i,2),RV(i,3),miu);
 eh_record(i)=eh;
 if(ialpha~=1) disp(['distance=',num2str(timedistance(i,2))]); end;
 capFc=acosh((md/ah+1)/eh);
 capmhc=eh*sinh(capFc)-capFc;
%  Nsubstep=ceil((capmhc-capmh)/rec_dcapm);
 dcapm=(capmhc-capmh)/Nsubstep;
 nhyper=(miu/ah^3)^(1/2);
 dt=dcapm/nhyper;
 else
     disp "mode input error : 1=hyper 2=straight line";
 end
 for j=1:(Nsubstep+1)
    if(mode==1)
     %hyper
     capmht=capmh+(j-1)*dcapm;
     [RRt(j,1),RRt(j,2),RRt(j,3),RVt(j,1),RVt(j,2),RVt(j,3)]=el2xv(miu,ialpha,ah,eh,inch,capomh,omegah,capmht);
    elseif(mode==2)
     %straight line
     RRt(j,:)=RR(i,:)+RV(i,:)*(j-1)*dt;
    end
     z=dot(RRt(j,:),h(i,:))/norm(h(i,:));
     FN=miu*z/norm(RRt(j,:))^3;
     dit(j,1)=FN*beta/na*coswf/(1+e*cosf)*180/pi*dt;
 end
di(i,1)=2*sum(dit(1:Nsubstep,1));
if(ialpha~=1) disp(['di_this_time=',num2str(di(i,1))]); end;
end
sigma_di=sum(di);
disp(['sigma_i=',num2str(sigma_di)]);
if(mode==1) disp(['eh=',num2str(eh)]); end;
clear a e inc capom om M beta na cosf coswf sinwf;
clear ah eh inch capomh capmh omegah M beta nhyper cosf coswf sinwf;
clear z;
clear i j;