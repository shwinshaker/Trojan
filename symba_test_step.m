
ffname='symba_intro_test';
fname={'test_fictiousTrojan';'test_fictiousTrojan_step+';'test_fictiousTrojan_step++'};
tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{1},'/follow.out']);
time=zeros(length(tpel),length(fname));
tpa=zeros(length(tpel),length(fname));
tpe=zeros(length(tpel),length(fname));
tpi=zeros(length(tpel),length(fname));
for i=1:length(fname)
    tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{i},'/follow.out']);
    time(:,i)=tpel(:,1);
    tpa(:,i)=tpel(:,3);
    tpe(:,i)=tpel(:,4);
    tpi(:,i)=tpel(:,5);
end

figure;
set(gcf,'Position',[400,100,700,500],'color','w');

subplot(3,1,1);
plot(0,30,'w+');hold all;
plot(time(:,1),tpa(:,1),'k-');
plot(time(:,2),tpa(:,2),'m--');
plot(time(:,3),tpa(:,3),'r-')
hold off;

subplot(3,1,2);
plot(0,0,'w+');hold all;
plot(time(:,1),tpe(:,1),'k-');
plot(time(:,2),tpe(:,2),'m--');
plot(time(:,3),tpe(:,3),'r-')
hold off;

subplot(3,1,3);
plot(0,10,'w+');hold all;
plot(time(:,1),tpi(:,1),'k-');
plot(time(:,2),tpi(:,2),'m--');
plot(time(:,3),tpi(:,3),'r-')
hold off;