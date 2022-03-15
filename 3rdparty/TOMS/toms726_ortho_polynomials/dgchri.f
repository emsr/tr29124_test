c
c
      subroutine dgchri(n,iopt,nu0,numax,deps,da,db,dx,dy,dalpha,dbeta,
     *nu,ierr,ierrc,dnu,drhor,drhoi,droldr,droldi,ds,ds0,ds1,ds2)
c
c This is a double-precision version of the routine  gchri.
c
      double precision deps,da(numax),db(numax),dx,dy,dalpha(n),
     *dbeta(n),dnu(*),drhor(*),drhoi(*),droldr(*),droldi(*),
     *ds(n),ds0(*),ds1(*),ds2(*)
c
c The arrays  dnu,drhor,drhoi,droldr,droldi,ds0,ds1,ds2  are assumed
c to have dimension  2*n.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      nd=2*n
      ndm1=nd-1
      if(iopt.eq.1) then
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 10 k=1,nd
          dnu(k)=-drhor(k)
   10   continue
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else if(iopt.eq.2) then
        dy=dabs(dy)
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 20 k=1,nd
          dnu(k)=-drhoi(k)/dy
   20   continue
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else
        ierr=1
        return
      end if
      end

