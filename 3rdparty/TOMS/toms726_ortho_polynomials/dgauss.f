c
c
      subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)
c
c This is a double-precision version of the routine  gauss.
c
      double precision dalpha,dbeta,deps,dzero,dweigh,de,dp,dg,dr,
     *ds,dc,df,db
      dimension dalpha(n),dbeta(n),dzero(n),dweigh(n),de(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      dzero(1)=dalpha(1)
      if(dbeta(1).lt.0.d0) then
        ierr=-2
        return
      end if
      dweigh(1)=dbeta(1)
      if (n.eq.1) return
      dweigh(1)=1.d0
      de(n)=0.d0
      do 100 k=2,n
        dzero(k)=dalpha(k)
        if(dbeta(k).lt.0.d0) then
          ierr=-2
          return
        end if
        de(k-1)=dsqrt(dbeta(k))
        dweigh(k)=0.d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(dabs(de(m)).le.deps*(dabs(dzero(m))+dabs(dzero(m+1)))) 
     *      goto 120
  110   continue
  120   dp=dzero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(dzero(l+1)-dp)/(2.d0*de(l))
        dr=dsqrt(dg*dg+1.d0)
        dg=dzero(m)-dp+de(l)/(dg+dsign(dr,dg))
        ds=1.d0
        dc=1.d0
        dp=0.d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*de(i)
          db=dc*de(i)
          if(dabs(df).lt.dabs(dg)) goto 150
          dc=dg/df
          dr=dsqrt(dc*dc+1.d0)
          de(i+1)=df*dr
          ds=1.d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=dsqrt(ds*ds+1.d0)
          de(i+1)=dg*dr
          dc=1.d0/dr
          ds=ds*dc
  160     dg=dzero(i+1)-dp
          dr=(dzero(i)-dg)*ds+2.d0*dc*db
          dp=ds*dr
          dzero(i+1)=dg+dp
          dg=dc*dr-db
          df=dweigh(i+1)
          dweigh(i+1)=ds*dweigh(i)+dc*df
          dweigh(i)=dc*dweigh(i)-ds*df
  200   continue
        dzero(l)=dzero(l)-dp
        de(l)=dg
        de(m)=0.d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=dzero(i)
        do 260 j=ii,n
          if(dzero(j).ge.dp) goto 260
          k=j
          dp=dzero(j)
  260   continue
        if(k.eq.i) goto 300
        dzero(k)=dzero(i)
        dzero(i)=dp
        dp=dweigh(i)
        dweigh(i)=dweigh(k)
        dweigh(k)=dp
  300 continue
      do 310 k=1,n
        dweigh(k)=dbeta(1)*dweigh(k)*dweigh(k)
  310 continue
      return
  400 ierr=l
      return
      end

