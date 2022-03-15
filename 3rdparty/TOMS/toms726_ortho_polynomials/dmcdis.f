c
c
      subroutine dmcdis(n,ncapm,mc,mp,dxp,dyp,dquad,deps,iq,idelta,
     *irout,finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncap,
     *kount,ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
c
c This is a double-precision version of the routine  mcdis.
c
      double precision dxp,dyp,deps,dendl,dendr,dxfer,dwfer,dalpha,
     *dbeta,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2
      dimension dxp(*),dyp(*),dendl(mc),dendr(mc),dxfer(ncapm),
     *dwfer(ncapm),dalpha(n),dbeta(n),dbe(n),dx(ncapm),dw(ncapm),
     *dxm(*),dwm(*),dp0(*),dp1(*),dp2(*)
      logical finld,finrd
c
c The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,
c the arrays  dxm,dwm,dp0,dp1,dp2  dimension  mc*ncapm+mp.
c
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierrd=-1
        return
      end if
      incap=1
      kount=-1
      ierr=0
      do 10 k=1,n
        dbeta(k)=0.d0
   10 continue
      ncap=(2*n-1)/idelta
   20 do 30 k=1,n
        dbe(k)=dbeta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierrd=ncapm
        return
      end if
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call dquad(ncap,dx,dw,i,ierr)
        else
          call dqgp(ncap,dx,dw,i,ierr,mc,finld,finrd,dendl,dendr,dxfer,
     *      dwfer)
        end if
        if(ierr.ne.0) then
          ierrd=i
          return
        end if
        do 40 k=1,ncap
          dxm(im1tn+k)=dx(k)
          dwm(im1tn+k)=dw(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          dxm(mtncap+k)=dxp(k)
          dwm(mtncap+k)=dyp(k)
   60   continue
      end if
      if(irout.eq.1) then
        call dsti(n,mtncap+mp,dxm,dwm,dalpha,dbeta,ied,dp0,dp1,dp2)
      else
        call dlancz(n,mtncap+mp,dxm,dwm,dalpha,dbeta,ied,dp0,dp1)
      end if
      do 70 k=1,n
        if(dabs(dbeta(k)-dbe(k)).gt.deps*dabs(dbeta(k))) goto 20
   70 continue
      return
      end

