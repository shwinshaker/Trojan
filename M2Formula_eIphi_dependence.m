% CE times - total inc
if exist('el','var') && exist('RCE','var') ...,
      && exist('stdCE','var') && exist('CE_times','var')
    tag='repeat';
else
    tag='data';
end
if ~strcmp(tag,'repeat')
clear;
Dir='ServerMount';

Dirname={'pl'};

filename=['~/Documents/',Dir,'/LAB/tools/Plutino/PlutinoList.txt'];
fid=fopen(filename,'r');
Pluinfo=textscan(fid,'%s %f %f %f %f %f %f %f %f %f');
fclose(fid);
PluPhi=Pluinfo{10};
pl_name=unique(Pluinfo{:,1});

Ndata=0;
for ii=1:length(Dirname)
    filename=eval([Dirname{ii},'_name']);
    Ndata=Ndata+length(filename);
end
CE_times=zeros(Ndata,1);
RCE=zeros(Ndata,1);
stdCE=zeros(Ndata,1);
el=zeros(Ndata,3);
name=cell(Ndata,1);

fontsize=15;
errid=[];
for ii=1:length(Dirname)
    disp(Dirname{ii});
    filename=eval([Dirname{ii},'_name']);
    if strcmp(Dirname{ii},'pl')
        ffname='RealPlutinos';
        di_name='di_record_perturb';
        CE_name='CE_record';
    elseif strcmp(Dirname{ii},'ran_npl')
        ffname='RanPlutinos';
        di_name='di_fit_perturb';
        CE_name='ran_record';
    elseif intersect(Dirname{ii},['npl' 'npl20' 'npl40'])
        ffname='RealPlutinosNpl';
        di_name='di_record_perturb';
        CE_name='CE_record';
    end
    for i=1:length(filename)
        
        if strcmp(Dirname{ii},'ran_npl')
            fname=filename{i};
        elseif strcmp(Dirname{ii},'npl20')
            fname=[filename{i},'_1Gyr_20pl'];
        elseif strcmp(Dirname{ii},'npl40')
            fname=[filename{i},'_1Gyr_40pl'];
        else
            fname=[filename{i},'_1Gyr'];
        end
        disp(filename{i});
        try
            plel=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/plel.txt'));
        catch
            errid=[errid;i];
            continue;
        end
        el(i,1)=mean(plel(:,2));
        if el(i,1)>40 || el(i,1)<39
            continue;
        end
        name{i}=filename{i};        
        di_record_perturb=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/',di_name,'.txt'));
        RCE(i)=sum((di_record_perturb-mean(di_record_perturb)).^2);
        stdCE(i)=(var(di_record_perturb))^(1/2);
        CE_times(i)=length(di_record_perturb);
        Nepel=load(strcat('~/Documents/',Dir,'/LAB/CE_realp/',ffname,'/',fname,'/Nepel.txt'));

        el(i,2)=mean(plel(:,3));
        el(i,3)=mean(plel(:,4));
        phi=mod(3*(plel(:,5)+plel(:,6)+plel(:,7))-2*(Nepel(:,5)+Nepel(:,6)+Nepel(:,7))-(plel(:,5)+plel(:,6)),360);
        %phi=mod(3*(Nepel(:,5)+Nepel(:,6)+Nepel(:,7))-2*(plel(:,5)+plel(:,6)+plel(:,7))-(Nepel(:,5)+Nepel(:,6)),360);
        %plela=load(strcat('~/Documents/',Dir,'/LAB/tools/Plutino/',filename{i},'.txt'));
        %phia=plela(9);
        lineid=find(strcmp(Pluinfo{:,1},filename{i}));
        if length(lineid)>1
            lineid=lineid(1);
        end
        phia=PluPhi(lineid);
        disp(mean(phi));
        el(i,4)=(max(phi)-min(phi))/2;
        el(i,5)=phia;
    end
end

%name{errid}=[];
el(cellfun(@isempty,name),:)=[];
CE_times(cellfun(@isempty,name))=[];
RCE(cellfun(@isempty,name),:)=[];
stdCE(cellfun(@isempty,name),:)=[];
name(cellfun(@isempty,name))=[];
namedata=[name num2cell(el) num2cell(CE_times)]; %% converted to cell

end

% figure(1);
% set(gcf,'Position',[400,100,1000,400],'color','w');
% 
% subplot(1,3,1);
% abnormalIndex=find(el(:,1)>40);
% plot(el(:,2),CE_times,'k.');hold on;
% plot(el(abnormalIndex,2),CE_times(abnormalIndex),'rx');
% subplot(1,3,2);
% plot(el(:,3),CE_times,'k.');hold on;
% plot(el(abnormalIndex,3),CE_times(abnormalIndex),'rx');
% subplot(1,3,3);
% plot(el(:,5),CE_times,'k.');hold on;
% plot(el(abnormalIndex,5),CE_times(abnormalIndex),'rx');

% figure(1);
% set(gcf,'Position',[400,100,1000,400],'color','w');
% 
% subplot(1,2,1);
% abnormalIndex=find(el(:,1)>40);
% plot(0,0,'w');hold on;
% color=[1-(CE_times+1)/max((CE_times+1)) 1-(CE_times+1)/max((CE_times+1)) 1-(CE_times+1)/max((CE_times+1))];
% scatter(el(:,2),el(:,3),50,color,'filled');
% %plot(el(abnormalIndex,2),CE_times(abnormalIndex),'rx');
% subplot(1,2,2);
% plot(0,0,'w');hold on;
% scatter(el(:,5),el(:,3),50,color,'filled');hold on;
% %plot(el(abnormalIndex,3),CE_times(abnormalIndex),'rx');

TNOFilename='~/Desktop/TNOData.txt';
disp(['TNO info file: ',TNOFilename]);
fid=fopen(TNOFilename,'r');
%headerlines=textscan(fid,'%s',1,'delimiter','\n');
info=textscan(fid,'%s %f %f %f %f %f %f');
fclose(fid);

TNOdata=[info{2:end}];%%% TNO Data
PlutinoData=TNOdata(TNOdata(:,1)>39 & TNOdata(:,1)<40,:);

figure(1);
set(gcf,'Position',[400,100,1000,400],'color','w');

subplot(1,3,1);
plot(0,0,'w');hold all;
for i=1:length(el)
    if CE_times(i)==0 
        color='r';
        marker='.';
        markersize=5;
%     elseif CE_times(i)<1000 && CE_times(i)>100
%         color='g';
%     elseif CE_times(i)>1000
%         color='r';
    else
        color='b';
        marker='o';
        markersize=CE_times(i)/max(CE_times)*20;
      %color=[log(CE_times(i))/log(max(CE_times)) 0 1-log(CE_times(i))/log(max(CE_times))];     
    end
    %plot(el(i,2),el(i,3),'.','color',color,'markersize',15);hold on;
    plot(el(i,2),el(i,3),[color,marker],'MarkerSize',markersize);hold on;grid on;

    %plot(PlutinoData(:,2),PlutinoData(:,3),'r.','markersize',5);
end
yylim=get(gca,'ylim');
sep=min(el(CE_times>0,2));
plot([sep sep],[yylim(1),yylim(2)],'k--');
set(text(sep+0.01,9/10*yylim(2),['$$e= ',num2str(sep),'$$']),'Interpreter','latex','fontsize',fontsize,'color','r');
hold off;
xlabel('$e$','fontsize',fontsize,'Interpreter','latex');
ylabel('$Inc.$','fontsize',fontsize,'Interpreter','latex');

subplot(1,3,2);
plot(0,0,'w');hold on;
% for i=1:length(el)
%     if CE_times(i)<100
%         color='b';
%     elseif CE_times(i)<1000 && CE_times(i)>100
%         color='g';
%     elseif CE_times(i)>1000
%         color='r';
%     end
%     plot(el(i,5),el(i,3),[color,'.'],'markersize',15);
% end
for i=1:length(el)
    if CE_times(i)==0 
        color='r';
        marker='.';
        markersize=5;
%     elseif CE_times(i)<1000 && CE_times(i)>100
%         color='g';
%     elseif CE_times(i)>1000
%         color='r';
    else
        color='b';
        marker='o';
        markersize=CE_times(i)/max(CE_times)*20;
      %color=[log(CE_times(i))/log(max(CE_times)) 0 1-log(CE_times(i))/log(max(CE_times))];     
    end
    %plot(el(i,2),el(i,3),'.','color',color,'markersize',15);hold on;
    plot(el(i,4),el(i,3),[color,marker],'MarkerSize',markersize);hold on;grid on;

    %plot(PlutinoData(:,2),PlutinoData(:,3),'r.','markersize',5);
end

subplot(1,3,3);
plot(0,0,'w');hold on;
% for i=1:length(el)
%     if CE_times(i)<100
%         color='b';
%     elseif CE_times(i)<1000 && CE_times(i)>100
%         color='g';
%     elseif CE_times(i)>1000
%         color='r';
%     end
%     plot(el(i,5),el(i,3),[color,'.'],'markersize',15);
% end
for i=1:length(el)
    if CE_times(i)==0 
        color='r';
        marker='.';
        markersize=5;
    else
        color='b';
        marker='o';
        markersize=CE_times(i)/max(CE_times)*20;
      %color=[log(CE_times(i))/log(max(CE_times)) 0 1-log(CE_times(i))/log(max(CE_times))];     
    end
    %plot(el(i,2),el(i,3),'.','color',color,'markersize',15);hold on;
    plot(el(i,4),el(i,2),[color,marker],'MarkerSize',markersize);hold on;grid on;

    %plot(PlutinoData(:,2),PlutinoData(:,3),'r.','markersize',5);
end
