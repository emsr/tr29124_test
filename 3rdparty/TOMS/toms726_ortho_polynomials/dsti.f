c
c
      subroutine dsti(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1,dp2)
c
c This is a double-precision version of the routine  sti.
c
      double precision dx,dw,dalpha,dbeta,dp0,dp1,dp2,dtiny,d1mach,
     *dhuge,dsum0,dsum1,dsum2,dt
      dimension dx(ncap),dw(ncap),dalpha(n),dbeta(n),dp0(ncap),
     *dp1(ncap),dp2(ncap)
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
      dsum0=0.d0
      dsum1=0.d0
      do 10 m=1,ncap
        dsum0=dsum0+dw(m)
        dsum1=dsum1+dw(m)*dx(m)
   10 continue
      dalpha(1)=dsum1/dsum0
      dbeta(1)=dsum0
      if(n.eq.1) return
      do 20 m=1,ncap
        dp1(m)=0.d0
        dp2(m)=1.d0
   20 continue
      do 40 k=1,nm1
        dsum1=0.d0
        dsum2=0.d0
        do 30 m=1,ncap
          if(dw(m).eq.0.d0) goto 30
          dp0(m)=dp1(m)
          dp1(m)=dp2(m)
          dp2(m)=(dx(m)-dalpha(k))*dp1(m)-dbeta(k)*dp0(m)
          if(dabs(dp2(m)).gt.dhuge .or. dabs(dsum2).gt.dhuge) then
            ierr=k
            return
          end if
          dt=dw(m)*dp2(m)*dp2(m)
          dsum1=dsum1+dt
          dsum2=dsum2+dt*dx(m)
   30   continue
        if(dabs(dsum1).lt.dtiny) then
          ierr=-k
          return
        end if
        dalpha(k+1)=dsum2/dsum1
        dbeta(k+1)=dsum1/dsum0
        dsum0=dsum1
   40 continue
      return
      end

