c
c
      subroutine dlancz(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1)
c
c This is a double-precision version of the routine  lancz.
c
      double precision dx(ncap),dw(ncap),dalpha(n),dbeta(n),
     *dp0(ncap),dp1(ncap),dpi,dgam,dsig,dt,dxlam,drho,dtmp,
     *dtsig,dtk
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        dp0(i)=dx(i)
        dp1(i)=0.d0
   10 continue
      dp1(1)=dw(1)
      do 30 i=1,ncap-1
        dpi=dw(i+1)
        dgam=1.d0
        dsig=0.d0
        dt=0.d0
        dxlam=dx(i+1)
        do 20 k=1,i+1
          drho=dp1(k)+dpi
          dtmp=dgam*drho
          dtsig=dsig
          if(drho.le.0.d0) then
            dgam=1.d0
            dsig=0.d0
          else
            dgam=dp1(k)/drho
            dsig=dpi/drho
          end if
          dtk=dsig*(dp0(k)-dxlam)-dgam*dt
          dp0(k)=dp0(k)-(dtk-dt)
          dt=dtk
          if(dsig.le.0.d0) then
            dpi=dtsig*dp1(k)
          else
            dpi=(dt**2)/dsig
          end if
          dtsig=dsig
          dp1(k)=dtmp
   20   continue
   30 continue
      do 40 k=1,n
        dalpha(k)=dp0(k)
        dbeta(k)=dp1(k)
   40 continue
      return
      end 

