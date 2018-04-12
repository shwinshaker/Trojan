clear;

ffname='RealPlutinosNpl';
fname='1999CE119_1Gyr_40pl';

ierror=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/ierror.txt']);

tag='di';

di_record_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt']);
di_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt']);

di_record_perturb(ierror)=[];
di_record_inout(ierror)=[];

di_fit_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt']);


tag='de';

de_record_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_perturb.txt']);
de_record_inout=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_record_inout.txt']);

de_record_perturb(ierror)=[];
de_record_inout(ierror)=[];

de_fit_perturb=load(['~/Documents/ServerMount/LAB/CE_realp/',ffname,'/',fname,'/',tag,'_fit_perturb.txt']);


err=(de_fit_perturb-de_record_inout)./de_record_inout;
%semilogx(de_record_inout,err,'k.');

Nbar=50;
maxerr=5;
errxx=1:Nbar;
errxx=errxx';
errxx=errxx/Nbar*maxerr;
errxx=[-errxx;0;errxx];
errxx=sort(errxx);

plotx=(errxx(2:end)+errxx(1:end-1))/2;
countx=histcounts(err,errxx)'/length(err);

plot(plotx,countx,'r-');


