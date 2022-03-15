c
c
      subroutine dmcheb(n,ncapm,mcd,mp,dxp,dyp,dquad,deps,iq,
     *idelta,finld,finrd,dendl,dendr,dxfer,dwfer,da,db,dnu,dalpha,
     *dbeta,ncap,kount,ierrd,dbe,dx,dw,dxm,dwm,ds,ds0,ds1,ds2)
c
c This is a double-precision version of the routine  mccheb.
c
      double precision dxp,dyp,deps,dendl,dendr,dxfer,dwfer,da,db,
     *dnu,dalpha,dbeta,dbe,dx,dw,dxm,dwm,ds,ds0,ds1,ds2,dsum,dp1,
     *dp,dpm1
      dimension dxp(*),dyp(*),dendl(mcd),dendr(mcd),dxfer(ncapm),
     *dwfer(ncapm),da(*),db(*),dnu(*),dalpha(n),dbeta(n),dbe(n),
     *dx(ncapm),dw(ncapm),dxm(*),dwm(*),ds(n),ds0(*),ds1(*),ds2(*)
      logical finld,finrd
c
c The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,
c the arrays  da,db  dimension 2*n-1, the arrays  dnu,ds0,ds1,ds2
c dimension  2*n, and the arrays  dxm,dwm  dimension  mc*ncapm+mp.
c
      nd=2*n
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierrd=-1
        return
      end if
      incap=1
      kount=-1
      ierrd=0
      do 10 k=1,n
        dbeta(k)=0.d0
   10 continue
      ncap=(nd-1)/idelta
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
      mtncap=mcd*ncap
      do 50 i=1,mcd
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call dquad(ncap,dx,dw,i,ierr)
        else
          call dqgp(ncap,dx,dw,i,ierr,mcd,finld,finrd,dendl,dendr,
     *      dxfer,dwfer)
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
      mtnpmp=mtncap+mp
      do 90 k=1,nd
        km1=k-1
        dsum=0.d0
        do 80 i=1,mtnpmp
          dp1=0.d0
          dp=1.d0
          if(k.gt.1) then
            do 70 l=1,km1
              dpm1=dp1
              dp1=dp
              dp=(dxm(i)-da(l))*dp1-db(l)*dpm1
   70       continue
          end if
          dsum=dsum+dwm(i)*dp
   80   continue
        dnu(k)=dsum
   90 continue
      call dcheb(n,da,db,dnu,dalpha,dbeta,ds,ierr,ds0,ds1,ds2)
      do 100 k=1,n
        if(dabs(dbeta(k)-dbe(k)).gt.deps*dabs(dbeta(k))) goto 20
  100 continue
      return
      end

