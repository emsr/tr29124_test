c
c
      subroutine dqgp(n,dx,dw,i,ierr,mcd,finld,finrd,dendl,dendr,
     *  dxfer,dwfer)
c
c This is a double-precision version of the routine  qgp. The user
c has to supply the routine
c
c              double precision function dwf(dx,i),
c
c which evaluates the weight function in double precision at the
c point  dx  on the i-th component interval.
c
      double precision dx,dw,dendl,dendr,dxfer,dwfer,dphi,dphi1,dwf
      dimension dx(n),dw(n),dendl(mcd),dendr(mcd),dxfer(*),dwfer(*)
      logical finld,finrd
c
c The arrays  dxfer,dwfer  are dimensioned in the routine  dmcdis.
c
      ierr=0
      if(i.eq.1) call dfejer(n,dxfer,dwfer)
      if(i.gt.1 .and. i.lt.mcd) goto 60
      if(mcd.eq.1) then
        if(finld.and.finrd) goto 60
        if(finld) goto 20
        if(finrd) goto 40
        do 10 k=1,n
          call dsymtr(dxfer(k),dphi,dphi1)
          dx(k)=dphi
          dw(k)=dwfer(k)*dwf(dphi,i)*dphi1
   10   continue
        return
      else
        if((i.eq.1.and.finld).or.(i.eq.mcd.and.finrd)) goto 60
        if(i.eq.1) goto 40
      end if
   20 do 30 k=1,n
        call dtr(dxfer(k),dphi,dphi1)
        dx(k)=dendl(mcd)+dphi
        dw(k)=dwfer(k)*dwf(dx(k),mcd)*dphi1
   30 continue
      return
   40 do 50 k=1,n
        call dtr(-dxfer(k),dphi,dphi1)
        dx(k)=dendr(1)-dphi
        dw(k)=dwfer(k)*dwf(dx(k),1)*dphi1
   50 continue
      return
   60 do 70 k=1,n
        dx(k)=.5d0*((dendr(i)-dendl(i))*dxfer(k)+dendr(i)+dendl(i))
        dw(k)=.5d0*(dendr(i)-dendl(i))*dwfer(k)*dwf(dx(k),i)
   70 continue
      return
      end

      subroutine dsymtr(dt,dphi,dphi1)
c
c This is a double-precision version of  symtr.
c
      double precision dt,dphi,dphi1,dt2
      dt2=dt*dt
      dphi=dt/(1.-dt2)
      dphi1=(dt2+1.d0)/(dt2-1.d0)**2
      return
      end

      subroutine dtr(dt,dphi,dphi1)
c
c This is a double-precision version of  tr.
c
      double precision dt,dphi,dphi1
      dphi=(1.d0+dt)/(1.d0-dt)
      dphi1=2.d0/(dt-1.d0)**2
      return
      end

      subroutine dfejer(n,dx,dw)
c
c This is a double-precision version of  fejer.
c
      double precision dx,dw,dpi,dn,dc1,dc0,dt,dsum,dc2
      dimension dx(n),dw(n)
      dpi=4.d0*datan(1.d0)
      nh=n/2
      np1h=(n+1)/2
      dn=dble(n)
      do 10 k=1,nh
        dx(n+1-k)=dcos(.5d0*dble(2*k-1)*dpi/dn)
        dx(k)=-dx(n+1-k)
   10 continue
      if(2*nh.ne.n) dx(np1h)=0.d0
      do 30 k=1,np1h
        dc1=1.d0
        dc0=2.d0*dx(k)*dx(k)-1.d0
        dt=2.d0*dc0
        dsum=dc0/3.d0
        do 20 m=2,nh
          dc2=dc1
          dc1=dc0
          dc0=dt*dc1-dc2
          dsum=dsum+dc0/dble(4*m*m-1)
   20   continue
        dw(k)=2.d0*(1.d0-2.d0*dsum)/dn
        dw(n+1-k)=dw(k)
   30 continue
      return
      end

