function [ialpha,a,e,inc,capom,omega,capm]=xv2el(x,y,z,vx,vy,vz,gmsum)
TINY=4e-15;
hx=y*vz-z*vy;
hy=z*vx-x*vz;
hz=x*vy-y*vx;
h2=hx*hx+hy*hy+hz*hz;
h=sqrt(h2);
inc=acos(hz/h);

%Compute longitude of ascending node CAPOM and the argument of latitude u.
fac = sqrt(hx^2 + hy^2)/h;
if(fac<TINY )
    capom=0;
	 u=atan2(y,x);
     if(abs(inc-pi)<10*TINY) 
          u=-u;  
     end
else
	  capom=atan2(hx,-hy);	  
	  u = atan2( z/sin(inc) , x*cos(capom) + y*sin(capom));
end
if(capom<0) 
    capom = capom + 2*pi; 
end
if(u<0) 
    u = u + 2*pi; 
end

r = sqrt(x*x + y*y + z*z);
v2 = vx*vx + vy*vy + vz*vz;
v = sqrt(v2);
vdotr = x*vx + y*vy + z*vz;
energy = 0.5*v2 - gmsum/r;

if(abs(energy*r/gmsum)<sqrt(TINY))
	   ialpha = 0;
else
    if(energy<0)
        ialpha = -1;  
    end
    if(energy>0)
        ialpha = +1; 
    end
end
%ELLIPSE :
if(ialpha==-1)
      disp('Warning: Ellipsic orbits detected!!!');
	  a = -0.5*gmsum/energy;  
	  fac = 1 - h2/(gmsum*a);
      if (fac>TINY)
             e = sqrt ( fac );
             face =(a-r)/(a*e);
%Apr. 16/93 : watch for case where face is slightly outside unity
             if ( face > 1 )  
                cape = 0 ;
             else
                 if ( face > -1 )
                   cape = acos( face );
                else
                   cape = pi;
                end
             end
             if ( vdotr < 0  )
                cape = 2 *pi - cape; 
             end
	    cw = (cos( cape) -e)/(1  - e*cos(cape));
	    sw = sqrt(1  - e*e)*sin(cape)/(1  - e*cos(cape));
	    w = atan2(sw,cw);
        if(w < 0 ) 
            w = w + 2 *pi;
        end
      else
	    e = 0; 
	    w = u;
	    cape = u;
      end
      capm = cape - e*sin (cape);
	  omega = u - w;
      if(omega < 0 ) 
          omega = omega + 2 *pi;
      end
	  omega = mod(omega,2*pi);  
end
      
 %HYPERBOLA
 if(ialpha == +1)  
        a = +0.5*gmsum/energy;  
        fac = h2/(gmsum*a);
      if (fac > TINY)  
 	    e = sqrt ( 1  + fac );
	    tmpf = (a+r)/(a*e);
        if(tmpf<1.0)  
            tmpf = 1.0;
        end
	    capf = log(tmpf + sqrt(tmpf*tmpf -1 ));
        if ( vdotr < 0  ) 
            capf = - capf; 
        end
	    cw = (e - cosh(capf))/(e*cosh(capf) - 1  );
	    sw = sqrt(e*e - 1 )*sinh(capf)/(e*cosh(capf) - 1  );
	    w = atan2(sw,cw);
        if(w < 0 ) 
            w = w + 2 *pi; 
        end
	  else
%we only get here if a hyperbola is essentially a parabola
%so we calculate e and w accordingly to avoid singularities
	    e = 1 ;
	    tmpf = 0.5*h2/gmsum;
	    w = acos(2 *tmpf/r -1 );
        if ( vdotr < 0 ) 
            w = 2 *pi - w; 
        end
	    tmpf = (a+r)/(a*e);
	    capf = log(tmpf + sqrt(tmpf*tmpf -1 ));
      end

	  capm = e * sinh(capf) - capf;
	  omega = u - w;
      if(omega < 0 ) 
          omega = omega + 2 *pi; 
      end
	  omega = mod(omega,2*pi); 
	end


% PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)
if(ialpha == 0)  
    disp('Warning: Parabolic orbits detected!!!');
	  a =  0.5*h2/gmsum  ;
	  e = 1 ;
	  w = acos(2 *a/r -1 );
      if ( vdotr < 0 ) 
          w = 2 *pi - w;
      end
	  tmpf = tan(0.5 * w);
	  capm = tmpf* (1  + tmpf*tmpf/3 );
	  omega = u - w;
      if(omega < 0 ) 
          omega = omega + 2 *pi;
      end
	  omega = mod(omega,2*pi); 
end
 %orbel_xv2el



