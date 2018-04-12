clear;
ffname='RealPlutinosNpl';
%ffname='RanPlutinos';
fname='1999CE119_1Gyr_40pl';
%fname='2001KN77_1Gyr_40pl';
%fname='NewRanCETest';

PV_record_pl=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_pl.txt'));
PV_record_tp=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/PV_record_tp.txt'));
di_record_inout=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/di_record_perturb.txt'));


r2hill_record=load(strcat('~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));
r2hill_mean=mean(r2hill_record);
hill=sqrt(r2hill_mean);

PV_rlt=PV_record_tp-PV_record_pl;
x=(PV_rlt(:,1).^2+PV_rlt(:,2).^2).^(1/2);
y=-PV_rlt(:,3);
theta=atan(y./x);

dis=(x.^2+y.^2).^(1/2);
maxdis=max(dis);
% 
figure(2);
plot(0,0,'w');hold all;
for i=1:length(theta)
    plot(theta(i),di_record_inout(i),'.','color',[dis(i)/maxdis 0 1-dis(i)/maxdis]);
end
hold off;
