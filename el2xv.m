function [x,y,z,vx,vy,vz]=el2xv(gm,ialpha,a,e,inc,capom,omega,capm)
TINY=4e-15;

% Executable code 
if(e<0)  
    disp( ' ERROR in orbel_el2xv: e<0, setting e=0!!1');
    e = 0;
end

% check for inconsistencies between ialpha and e
em1 = e - 1;
if( ((ialpha ==0) && (abs(em1)>TINY))  ||((ialpha<0) && (e>1.0 ))  ||((ialpha>0) && (e<1.0 )) )   
            disp( 'ERROR in orbel_el2xv: ialpha and e inconsistent');
            disp( '                       ialpha = ',ialpha);
            disp( '                            e = ',e);
end

%Generate rotation matrices (on p. 42 of Fitzpatrick)

	sp=sin(omega);cp=cos(omega);
	so=sin(capom);co=cos(capom);
	si=sin(inc);ci=cos(inc);
	d11 = cp*co - sp*so*ci;
	d12 = cp*so + sp*co*ci;
	d13 = sp*si;
	d21 = -sp*co - cp*so*ci;
	d22 = -sp*so + cp*co*ci;
	d23 = cp*si;
    
%Get the other quantities depending on orbit type ( i.e. IALPHA)
if (ialpha  == -1)  
	  cape = ehybrid(e,capm);
	  scap=sin(cape);ccap=cos(cape);
	  sqe = sqrt(1  -e*e);
	  sqgma = sqrt(gm*a);
	  xfac1 = a*(ccap - e);
	  xfac2 = a*sqe*scap;
	  ri = 1 /(a*(1  - e*ccap));
	  vfac1 = -ri * sqgma * scap;
	  vfac2 = ri * sqgma * sqe * ccap;
end

if (ialpha  == +1)
	  capf = fhybrid(e,capm);
	  shcap=sinh(capf);chcap=cosh(capf);
	  sqe = sqrt(e*e - 1  );
	  sqgma = sqrt(gm*a);
	  xfac1 = a*(e - chcap);
	  xfac2 = a*sqe*shcap;
	  ri = 1 /(a*(e*chcap - 1 ));
	  vfac1 = -ri * sqgma * shcap;
	  vfac2 = ri * sqgma * sqe * chcap;
end
if (ialpha  == 0)
	  zpara = zget(capm);
	  sqgma = sqrt(2 *gm*a);
	  xfac1 = a*(1  - zpara*zpara);
	  xfac2 = 2 *a*zpara;
	  ri = 1 /(a*(1  + zpara*zpara));
	  vfac1 = -ri * sqgma * zpara;
	  vfac2 = ri * sqgma ;
end
    
	x =  d11*xfac1 + d21*xfac2;
	y =  d12*xfac1 + d22*xfac2;
	z =  d13*xfac1 + d23*xfac2;
	vx = d11*vfac1 + d21*vfac2;
	vy = d12*vfac1 + d22*vfac2;
	vz = d13*vfac1 + d23*vfac2;
        % orbel_el2xv
function E=ehybrid(e,m)
if(e < 0.18d0) 
	  E = esolmd(e,m);
else
    if( e <= 0.8d0)
	     E = eget(e,m);
    else
	     E = ehie(e,m);
    end
end
return

function F=fhybrid(e,n)
abn = n;
if(n<0) 
    abn = -abn;
end
if(abn < 0.636*e-0.6) 
	  F = flon(e,n);
else 
	  F = fget(e,n);
end
return

function Z=zget(q)
iflag = 0;
if(q<0)  
    iflag = 1;
    q = -q;
end
if (q<1e-3)  
    Z = q*(1.d0 - (q*q/3.d0)*(1 -q*q));
else
    x = 0.5*(3*q + sqrt(9*(q^2) +4));
    tmp = x^(1/3);
    Z = tmp - 1/tmp;
end
if(iflag ==1)  
    Z = -Z;
    q = -q;
end
return

function eget=eget(e,m)
sm=sin(m);cm=cos(m);
% c  begin with a guess accurate to order ecc**3	
x = m + e*sm*( 1 + e*( cm + e*( 1 -1.5*sm*sm)));
% c  Go through one iteration for improved estimate
sx=sin(x);cx=cos(x);
	  es = e*sx;
	  ec = e*cx;
	  f = x - es  - m;
	  fp = 1.d0 - ec ;
	  fpp = es ;
	  fppp = ec ;
	  dx = -f/fp;
	  dx = -f/(fp + dx*fpp/2.d0);
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0);
      eget = x + dx;
% c Do another iteration.
% c For m between 0 and 2*pi this seems to be enough to
% c get near roundoff error for eccentricities between 0 and 0.8
 x =eget;
 sx=sin(x);cx=cos(x);
 es = e*sx;
 ec = e*cx;
 f = x - es  - m;
 fp = 1 - ec ;
 fpp = es ;
 fppp = ec ;
 dx = -f/fp;
 dx = -f/(fp + dx*fpp/2.d0);
 dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0);
eget = x + dx;
return

function esolmd=esolmd(e,m)
	 sm=sin(m);cm=cos(m);
	  x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)));
	  sx=sin(x);cx=cos(x);
	  es = e*sx;
	  ec = e*cx;
	  f = x - es  - m;
	  fp = 1.d0 - ec ;
	  fpp = es ;
	  fppp = ec ;
	  dx = -f/fp;
	  dx = -f/(fp + dx*fpp/2.d0);
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0);
esolmd = x + dx;
return

function ehie=ehie(e,m)
NMAX = 3;
TWOPI=2*pi;
	iflag = 0;
	nper = m/TWOPI;
	m = m - nper*TWOPI;
	if (m<0.d0) 
        m = m + TWOPI; 
    end

	if (m>pi) then
	   m = TWOPI - m;
	   iflag = 1;
    end
% 
% c Make a first guess that works well for e near 1.
	x = (6.d0*m)^(1.d0/3.d0) - m;
	niter =0;
% 
% c Iteration loop
for niter =1:NMAX
	   sa=sin(x + m);ca=cos(x+m);
	    esa = e*sa;
	    eca = e*ca;
	    f = x - esa;
	    fp = 1.d0 -eca;
	    dx = -f/fp;
	    dx = -f/(fp + 0.5d0*dx*esa);
	    dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx));
	    x = x + dx;
end
      ehie = m + x;
if (iflag==1) 
	  ehie = TWOPI - ehie;
	  m = TWOPI - m;
end
return

function fget=fget(e,capn)
IMAX = 10;
TINY=4e-15;
if( capn < 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0;
	   x = -log(tmp);
else
	   tmp = +2.d0*capn/e + 1.8d0;
	   x = log( tmp);
end
	fget = x;
for i = 1:IMAX
	  shx=sinh(x);chx=cosh(x);
	  esh = e*shx;
	  ech = e*chx;
	  f = esh - x - capn;
% c	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  ;
	  fpp = esh ;
	  fppp = ech ;
	  dx = -f/fp;
	  dx = -f/(fp + dx*fpp/2.d0);
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0);
	  fget = x + dx;
% c   If we have converged here there's no point in going on
	  if(abs(dx)<= TINY) 
          return
      end
	  x = fget;
end	
disp( 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' );
return

function flon=flon(e,capn)
TINY=4e-15;
IMAX = 10;
a11 = 156.d0;
a9 = 17160.d0;
a7 = 1235520.d0;
a5 = 51891840.d0;
a3 = 1037836800.d0;
b11 = 11.d0*a11;
b9 = 9.d0*a9;
b7 = 7.d0*a7;
b5 = 5.d0*a5;
b3 = 3.d0*a3;
	iflag = 0;
    if( capn < 0.d0) 
	   iflag = 1;
	   capn = -capn;
    end

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e);
	a0 = -6227020800.d0*capn/e;
	b1 = a1;
% 
% c  Set iflag nonzero if capn < 0., in which case solve for -capn
% c  and change the sign of the final answer for F.
% c  Begin with a reasonable guess based on solving the cubic for small F	
	a = 6.d0*(e-1.d0)/e;
	b = -6.d0*capn/e;
	sq = sqrt(0.25*b*b +a*a*a/27.d0);
	biga = (-0.5*b + sq)^0.3333333333333333d0;
	bigb = -(+0.5*b + sq)^0.3333333333333333d0;
	x = biga + bigb;
% c	write(6,*) 'cubic = ',x**3 +a*x +b
	flon = x;
% c If capn is tiny (or zero) no need to go further than cubic even for
% c e =1.
if( capn < TINY) 
        if(iflag == 1)
	   flon = -flon;
	   capn = -capn;
        end
        return
end

for i = 1:IMAX
	  x2 = x*x;
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))));
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))));
	  dx = -f/fp;
% c	  write(6,*) 'i,dx,x,f : '
% c	  write(6,432) i,dx,x,f
%432	  format(1x,i3,3(2x,1p1e22.15))
	 flon = x + dx;
% c   If we have converged here there's no point in going on
if(abs(dx)<=TINY) 
          if(iflag == 1)
              flon = -flon;
              capn = -capn;
          end
          return
end
	  x = flon;
end
% 
% c Abnormal return here - we've gone thru the loop 
% c IMAX times without convergence
if(iflag ==1) 
	   flon = -flon;
	   capn = -capn;
end
	disp( 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' )
	  diff = e*sinh(flon) - flon - capn;
	  disp( 'N, F, ecc*sinh(F) - F - N : ')
	  disp(capn,orbel_flon,diff)
return



