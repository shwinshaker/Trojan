%% ??????????????????? ??????
clear;

ffname='fit_record_contrast';
fname='1999CE119_1Gyr_5.0_hill';

CE_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/CE_record.txt'));
di_record_inout=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/','di_record_inout.txt'));
r2hill_record=load(strcat('~/Documents/swiftdata/LAB/CE_realp/',ffname,'/',fname,'/r2hill_record.txt'));

hill=sqrt(r2hill_record);
maxhill=5.0;
maxdistance=maxhill*hill;

distance=CE_record(:,2);
N=200;
Di=zeros(1,N);
for i=1:N
    Di(i)=sum(di_record_inout((i-1)/N*maxdistance<distance & distance<=i/N*maxdistance));
end

figure(1);
x=1:N;
x=x-1/2;
plot(x/N*maxhill,abs(Di),'k-');