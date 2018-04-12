%%need the closest position and velocity
%%need the r2hill info

clear;

fname='1999CE119_1Gyr_negative_revised_bin';
load(strcat('~/Documents/swiftdata/LAB/CE_realp/fit_record_contrast/',fname,'/ran_CE_di.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/fit_record_contrast/',fname,'/r2hill_record.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/fit_record_contrast/',fname,'/ran_CE_tp.txt'));
load(strcat('~/Documents/swiftdata/LAB/CE_realp/fit_record_contrast/',fname,'/ran_CE_pl.txt'));

vlim=2e-3;
md=3.5*r2hill_record^(1/2);

figure;
subplot(2,2,1);
theta=0:pi/360:2*pi;
x=md*cos(theta); 
y=md*sin(theta); 
plot(x,y,'r-');
axis square;axis([-md md -md md]);
ylabel('z /AU');
xlabel('(x^2+y^2)^{(1/2)} /AU');
hold on;
plot([0 0],[-md md],'r-');
hold on;
plot([-md md],[0 0],'r-');
hold on;
for i=1:length(ran_CE_tp)
plot(((ran_CE_tp(i,1)-ran_CE_pl(i,1))^2+(ran_CE_tp(i,2)-ran_CE_pl(i,2))^2)^(1/2),(ran_CE_tp(i,3)-ran_CE_pl(i,3)), ...,
'k.');
hold on;
end

subplot(2,2,2);
plot(0,0,'k+');
axis([-vlim vlim -vlim vlim]);axis square;hold on;
ylabel('Vz /AU/DAY');
xlabel('(Vx^2+Vy^2)^{(1/2)} /AU/DAY');
for i=1:length(ran_CE_tp)
plot(((ran_CE_tp(i,4)-ran_CE_pl(i,4))^2+(ran_CE_tp(i,5)-ran_CE_pl(i,5))^2)^(1/2),(ran_CE_tp(i,6)-ran_CE_pl(i,6)), ...,
'k.');
hold on;
end

subplot(2,2,4);
RVN=zeros(length(ran_CE_tp),1);
for i=1:length(ran_CE_tp)
    RVN(i,:)=norm([ran_CE_tp(i,4)-ran_CE_pl(i,4),ran_CE_tp(i,5)-ran_CE_pl(i,5),ran_CE_tp(i,6)-ran_CE_pl(i,6)]);
end
plot(ran_CE_di(:,1)/md*3.5,RVN(:),'k.');axis square;axis([0 3.5 0 vlim]);
ylabel('Vnorm /AU/DAY');
xlabel('CE distance /RHill');

clear i x y theta;
%%
%%%%3D
% [x,y,z]=sphere;mesh(md*x,md*y,md*z);
% alpha(0);
% axis([-0.18 0.18 -0.18 0.18 -0.18 0.18]);
% box on;axis square;grid on;grid minor;
% hold on;
% for i=1:length(ran_CE_tp)
% quiver3(ran_CE_tp(i,1)-ran_CE_pl(i,1),ran_CE_tp(i,2)-ran_CE_pl(i,2),ran_CE_tp(i,3)-ran_CE_pl(i,3),...,
%     10*(ran_CE_tp(i,4)-ran_CE_pl(i,4)),10*(ran_CE_tp(i,5)-ran_CE_pl(i,5)),10*(ran_CE_tp(i,6)-ran_CE_pl(i,6)),...,
%     'k','MaxHeadSize',1.0);
% hold on;
% pause(0.001); 
% end



