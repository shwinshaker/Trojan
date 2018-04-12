clear;
filename='~/Desktop/TNOData.txt';
disp(['TNO info file: ',filename]);
fid=fopen(filename,'r');
%headerlines=textscan(fid,'%s',1,'delimiter','\n');
info=textscan(fid,'%s %f %f %f %f %f %f');
fclose(fid);

data=[info{2:end}];

fontsize=15;

figure(1);
set(gcf,'Position',[400,100,500,800],'color','w');

subplot(2,1,1);
plot(40,0,'w');hold all;
markersize=(data(:,2)+1e-4)*100;
scatter(data(:,1),data(:,3),markersize,'k','filled');
hold off;
xlabel('a','fontsize',fontsize);
ylabel('Inc.','fontsize',fontsize)


subplot(2,1,2);
plot(40,0,'w');hold all;
markersize=(data(:,3)+1e-4);
scatter(data(:,1),data(:,2),markersize,'k','filled');
hold off;
xlabel('a','fontsize',fontsize);
ylabel('e','fontsize',fontsize)