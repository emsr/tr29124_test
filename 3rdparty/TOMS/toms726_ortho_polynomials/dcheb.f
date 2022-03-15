c
c
      subroutine dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
c
c This is a double-precision version of the routine  cheb.
c
      double precision da,db,dnu,dalpha,dbeta,ds,ds0,ds1,ds2,dtiny,
     *d1mach,dhuge
      dimension da(*),db(*),dnu(*),dalpha(n),dbeta(n),ds(n),
     *ds0(*),ds1(*),ds2(*)
c
c The arrays  da,db  are assumed to have dimension  2*n-1, the arrays
c dnu,ds0,ds1,ds2  dimension  2*n.
c
      nd=2*n
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      iderr=0
      if(dabs(dnu(1)).lt.dtiny) then
        iderr=1
        return
      end if
      if(n.lt.1) then
        iderr=2
        return
      end if
      dalpha(1)=da(1)+dnu(2)/dnu(1)
      dbeta(1)=dnu(1)
      if(n.eq.1) return
      ds(1)=dnu(1)
      do 10 l=1,nd
        ds0(l)=0.d0
        ds1(l)=dnu(l)
   10 continue
      do 40 k=2,n
        lk=nd-k+1
        do 20 l=k,lk
          ds2(l)=ds1(l+1)-(dalpha(k-1)-da(l))*ds1(l)-dbeta(k-1)*ds0(l)
     *      +db(l)*ds1(l-1)
        if(l.eq.k) ds(k)=ds2(k)
   20   continue
        if(dabs(ds(k)).lt.dtiny) then
          iderr=-(k-1)
          return
        else if(dabs(ds(k)).gt.dhuge) then
          iderr=k-1
          return
        end if
        dalpha(k)=da(k)+(ds2(k+1)/ds2(k))-(ds1(k)/ds1(k-1))
        dbeta(k)=ds2(k)/ds1(k-1)
        do 30 l=k,lk
          ds0(l)=ds1(l)
          ds1(l)=ds2(l)
   30   continue
   40 continue
      return
      end

