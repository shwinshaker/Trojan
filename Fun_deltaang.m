function [ deltaang ] = Fun_deltaang(aP,eP,IP,aT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% ????? ??????????? ????????
%% ?????????? ???????
%% ?????????????????? ???????? ?????= =?
% [~,y]=solve('(x/cosd(IP)+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
%     'x^2+y^2=aT^2','x','y');
[~,y]=solve('(x+aP*eP)^2+y^2/(1-eP^2)=aP^2',...,
    'x^2+y^2=aT^2','x','y');

yyy=0;
yy=eval(y);
for i=1:length(yy)
    if isreal(yy(i)) && yy(i) > 0
        yyy=yy(i);
    end
end
% disp(yyy);
angt=atand(yyy/(aT^2-yyy^2)^(1/2));
% disp(angt);
% angp=atand(yyy*cosd(IP)/((eP^2-1)*(aP^2*(eP^2-1)+yyy^2))^(1/2));
angp=atand(yyy/((eP^2-1)*(aP^2*(eP^2-1)+yyy^2))^(1/2));

% disp(angp);
deltaang=abs(angt-angp);

end

