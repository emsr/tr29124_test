c
c
      subroutine dkern(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,
     *  nu,ierr,droldr,droldi)
c
c This is a double-precision version of the routine  kern.
c
      double precision dx,dy,deps,da(numax),db(numax),dkerr(*),
     *  dkeri(*),droldr(*),droldi(*),dp0r,dp0i,dpr,dpi,dpm1r,
     *  dpm1i,dden,dt
c
c The arrays  dkerr,dkeri,droldr,droldi  are assumed to have
c dimension  n+1.
c
      call dknum(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,nu,ierr,
     *  droldr,droldi)
      if(ierr.ne.0) return
      dp0r=0.d0
      dp0i=0.d0
      dpr=1.d0
      dpi=0.d0
      do 10 k=1,n
        dpm1r=dp0r
        dpm1i=dp0i
        dp0r=dpr
        dp0i=dpi
        dpr=(dx-da(k))*dp0r-dy*dp0i-db(k)*dpm1r
        dpi=(dx-da(k))*dp0i+dy*dp0r-db(k)*dpm1i
        dden=dpr**2+dpi**2
        dt=(dkerr(k+1)*dpr+dkeri(k+1)*dpi)/dden
        dkeri(k+1)=(dkeri(k+1)*dpr-dkerr(k+1)*dpi)/dden
        dkerr(k+1)=dt
   10 continue
      return
      end

