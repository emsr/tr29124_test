c
c

      subroutine dlob(n,dalpha,dbeta,dleft,dright,dzero,dweigh,
     *ierr,de,da,db)
c
c This is a double-precision version of the routine  lob.
c
      double precision dleft,dright,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l,
     *dpm1r,ddet,dalpha(*),dbeta(*),dzero(*),dweigh(*),de(*),da(*),
     *db(*),d1mach
c
c The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
c dimension  n+2.
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0l=0.d0
      dp0r=0.d0
      dp1l=1.d0
      dp1r=1.d0
      do 20 k=1,np1
        dpm1l=dp0l
        dp0l=dp1l
        dpm1r=dp0r
        dp0r=dp1r
        dp1l=(dleft-da(k))*dp0l-db(k)*dpm1l
        dp1r=(dright-da(k))*dp0r-db(k)*dpm1r
   20 continue
      ddet=dp1l*dp0r-dp1r*dp0l
      da(np2)=(dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
      db(np2)=(dright-dleft)*dp1l*dp1r/ddet
      call dgauss(np2,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end

