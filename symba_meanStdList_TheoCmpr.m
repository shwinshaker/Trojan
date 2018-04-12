
clear;

ffname='symbaRealPlutinosNpl_fast';
fname={'1999CE119_2004UP10';'2001FU172_2004UP10';'1999CE119_2006RJ103';'2001FU172_2006RJ103'};
ffname1='symba_RealPlutinos';
fname1=fname;
titlename={'1999CE119&2004UP10';'2001FU172&2004UP10';'1999CE119&2006RJ103';'2001FU172&2006RJ103'};

% ffname1='RealPlutinos';
% fname1={'1999CE119_1Gyr';'2001FU172_1Gyr';'1999CE119&2006RJ103_1Gyr';'2001FU172&2006RJ103_1Gyr'};

% diname='de_record_inout';
% diname='di_record_inout';
diname='di_fit_perturb';
namesplit=strsplit(diname,'_');
prefix=namesplit{1};

inda=3;
Nbin=100;



if ~exist('meanStdList','var')
    
    meanStdList=zeros(length(fname),2,2);
    %% name * (mean or std) * (expr or theo)
    eiList=zeros(length(fname),2);
    %% name * (e or i)
    
    for isub=1:length(fname)
        
        di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname{isub},'/',diname,'.txt']);
        NCE=length(di_record_inout);
        
        di_norm=di_record_inout;%/Maxdi;
        di_norm_0=di_norm(di_norm~=0);
        diAbs=abs(di_norm);
        
        %% Theo
        plel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname1{isub},'/plel.txt']);
        tpel=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname1,'/',fname1{isub},'/tpel.txt']);
        
        pla=plel(:,inda);plamean=mean(pla);
        ple=plel(:,inda+1);plemean=mean(ple);%plemax=max(ple);plemin=min(ple);
        plinc=plel(:,inda+2);plincmean=mean(plinc);%plincstd=std(plinc);
        
        tpa=tpel(:,inda);tpamean=mean(tpa);
        tpe=tpel(:,inda+1);tpemean=mean(tpe);
        tpinc=tpel(:,inda+2);tpincmean=mean(tpinc);
        
        %N=length(di_record_inout); %% sample size
        %std0=abs(plincmean)^(1/3)+1;
        if strcmp(prefix,'di')
            [dinc,~,~]=Fun_diDstb_theo(1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
        else
            [~,dinc,~]=Fun_diDstb_theo(1,plamean,plemean,plincmean,tpamean,tpemean,tpincmean);
        end
        diTheoAbs=abs(dinc);
        
        meanStdList(isub,1,1)=mean(diAbs);
        meanStdList(isub,1,2)=mean(diTheoAbs);
        meanStdList(isub,2,1)=std(diAbs);
        meanStdList(isub,2,2)=std(diTheoAbs);
        
        H=kstest2(di_norm,dinc);
        disp('Theo-Numer test');
        disp(H);
        
    end
end

figure(1);
fontsize=15;
if strcmp(prefix,'di')
    xxlim=[1e-11 1];yylim=[1e-5 1];
else
    xxlim=[1e-11 0.1];yylim=[1e-5 1];
end
set(gcf,'Position',[400,100,700,500],'color','w');

plot(meanStdList(:,1,1),meanStdList(:,2,1),'k.');hold on;
plot(meanStdList(:,1,2),meanStdList(:,2,2),'r.');hold off;

