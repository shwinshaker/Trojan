
clear;
hillmax=3.5;
hill=0.0513;
dismax=hillmax*hill;
i=0;
Nrp=10;Nclm=20;maxpower=6;
% abszzsum=zeros(Nrp,Nclm);
% abstxsum=zeros(Nrp,Nclm);
% abstxmean=zeros(Nrp,Nclm);
abstxsum2=zeros(Nrp,Nclm);
abstxmean2=zeros(Nrp,Nclm);
txstd2=zeros(Nrp,Nclm);

for irp=1:Nrp
for i=1:Nclm
    iNran=1+maxpower/Nclm*(i-1);
    Nran=round(power(10,iNran));
%     dix=rand(Nran,1)*3.5*hill;
%     tx=rand(Nran,1)*2-1;
    tx2=randn(Nran,1);
%     zzran=(1-(dix/dismax).^2).^(1/2)./dix.*tx;
%     abszzsum(irp,i)=abs(sum(zzran));
%     abstxsum(irp,i)=abs(sum(tx));
%     abstxmean(irp,i)=abs(mean(tx));
    abstxsum2(irp,i)=abs(sum(tx2));
    abstxmean2(irp,i)=abs(mean(tx2));
    txstd2(irp,i)=var(tx2);
end
end
% abszzsummean=zeros(Nclm,1);
% abstxsummean=zeros(Nclm,1);
% abstxmeanmean=zeros(Nclm,1);
abstxsummean2=zeros(Nclm,1);
abstxmeanmean2=zeros(Nclm,1);
txstdmean2=zeros(Nclm,1);
for i=1:Nclm
%     abszzsummean(i)=trimmean(abszzsum(:,i),10);
%     abstxsummean(i)=trimmean(abstxsum(:,i),10);
%     abstxmeanmean(i)=trimmean(abstxmean(:,i),10);
    abstxsummean2(i)=trimmean(abstxsum2(:,i),10);
    abstxmeanmean2(i)=trimmean(abstxmean2(:,i),10);
    txstdmean2(i)=sum(txstd2(:,i));
end

%plotyy(abstxsummean,'k-',abstxmeanmean,'k--');hold on;
% plotyy(abstxmeanmean2,'k--',abstxsummean2,'k-');
plotyy(abstxsummean2,'k--',txstdmean2,'k-.');

% plot(txstdmean2,'k-.')
% plot(txsum,'b-')