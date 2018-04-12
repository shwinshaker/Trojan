function [countx,plotx,dx] = Fun_bin(xxlim1,xxlim2,Nbin)
        countx=(0:Nbin)';
        dx=(xxlim2-xxlim1)/Nbin;
        countx=xxlim1+countx*dx;
        plotx=(countx(2:end)-dx/2);
end
