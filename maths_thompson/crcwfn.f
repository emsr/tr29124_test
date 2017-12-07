cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine crcwfn(nch,ak2,kmax,coupl,rin,rout,
     *    accrcy,v,fcc,gcc,mch,maxv,correc,deriv,cfg,switch)
c
c     coupled real coulomb wavefunctions.
c     -----------------------------------
c     J A Christley, I J Thompson (Surrey) 1993
c
c     subroutine returns coupled coulomb wavefunctions (and derivs)
c     at rin integrated from boundary conditions at rout using the
c     potential defined by coupl.
c
c
      implicit real*8(a-h,o-z)
      real*8 ak2(mch),coupl(mch,mch,kmax),
     *    fcc(mch,mch,2),gcc(mch,mch,2),v(maxv),
     *    accrcy,rin(2),rout
      integer cstart,show
      logical correc,deriv,cfg
c
      show=0
      cstart=1
c
      do 10 i=1,maxv
 10       v(i)=0.0d0
c
      if (cfg) then
	al1= sqrt(coupl(1,1,2))+0.001d0
        lmax1=int(al1) + 30
        call coul(mch,nch,ak2,fcc,gcc,rout,v(1),
     *      v(lmax1+1),v(2*lmax1+1),v(3*lmax1+1),
     *      v(4*lmax1+1),v(4*lmax1+mch+1),lmax1,coupl,kmax)
        do 21  i=1,lmax1*4+5
21         v(i) = 0.0d0
      endif
c
      call cutup(nch,kmax,coupl,rin,rout,accrcy,v,deriv,mch,
     *           nsteps,nstepr2,maxv,switch,show)
c-----------------------------------------pointers for v()
      nch2 = nch*nch
      l1 = nsteps + 2
      l2 = l1 + 4*nch + 1
      imax = 2*nch*nsteps
      l3 = l2 + 4*nch2*(2*nsteps-1) +1
      l4 = l3 + nch+1
      l5 = l4 + nch+1
      l6 = l5 + nch2+1
      l7 = l6 + nch2+1
      l8 = l7 + nch2+1
      l9 = l8 + imax+1
      l10 = l9 + imax+1
      l11 = l10 + 12
      l12 = l11 + 2*nch +1
      if(show.ge.1) write(6,*)'no. of steps ',nsteps,nstepr2,deriv,imax
c
      if (l12.lt.maxv) then
c
      call matrix(nch,ak2,kmax,coupl,fcc,gcc,mch,
     *     correc,deriv,nsteps,nstepr2,v(1),v(l1),
     *     v(l2),l3-l2-1,v(l3),v(l4),v(l5),v(l6),v(l7),
     *     v(l8),v(l9),v(l10),imax,v(l11),switch,
     *     cstart,maxv,show)
c
      else
         write(6,*)'maxv array too small: need ',l12,' & have ',maxv
         stop
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cutup(nch,kmax,coupl,rin,rout,accrcy,v,
     *           deriv,mch,nsteps,nstepr2,maxv,switch,show)
c     devide integration range into piecewise segments
      implicit real*8(a-h,o-z)
      real*8 coupl(mch,mch,kmax),v(maxv),accrcy,rin(2),rout
      logical deriv,flagr2
      integer show
c
      flagr2 = .not. deriv
      if (abs(rin(1)-rin(2)).lt.1d-6)  then
         flagr2=.false.
         nstepr2=0
      endif
c
c---------------------------- divide integration range------------
      r = rin(1)
      ho = 0.1d0*r
      nsteps =0
      v(nsteps+1)=r
c
 10   nsteps = nsteps+1
      ddvm = 1.0d-30
      if (r.lt.switch) then
c        airy functions step depends on 2nd derivative
         do 20 ic=1,nch
            ddv = 0.0d0
            do 30 k=1,kmax
 30            ddv = ddv +  k*(k+1) * coupl(ic,ic,k)
     *           /(r**(k+2))
 20         if (abs(ddv).gt.ddvm) ddvm = abs(ddv)
c
         h = 4.0d0* sqrt(accrcy/ddvm)
      else
c        sines steps depend on 1st derivative
         do 40 ic=1,nch
            ddv = 0.0d0
            do 50 k=1,kmax
50             ddv = ddv - k * coupl(ic,ic,k)
     *           /(r**(k+1))
40          if (abs(ddv).gt.ddvm) ddvm = abs(ddv)
c        arbitrary factor 16
         h = 1.6d1*accrcy/ddvm
      endif
c
c     ensure next step is less than 1.4 * last
      if (h.gt.ho) h= ho
      ho = h*1.4d0
c
      if (flagr2) then
          if (r+h.ge.rin(2).and.rin(2).gt.rin(1)) then
              h =rin(2)-r
              nstepr2 = nsteps+1
              flagr2 = .false.
          endif
      endif
c
c   ensure at least two steps:
      if(nsteps.eq.1.and.r+h.gt.rout) h=(rout-rin(1))*0.5d0
c
      r = r+h
      if (r.lt.rout) then
c         ensure that we don't get a tiny last step
          if (rout-r.lt.0.6d0*h) r= rout -0.6d0*h
          v(nsteps+1) = r
          goto 10
      endif
      v(nsteps+1) = rout
c      if(show.ge.2) write(6,99) nsteps,rin,rout,(v(i),i=1,nsteps+1)
c99    format(' cutup: ',i5,' steps from',2f10.3,' to',f10.3,', at',
c     x      ,/,(1x,12f10.3))
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine matrix(nch,ak2,kmax,coupl,fcc,gcc,
     *      mch,correc,deriv,nsteps,nstepr2,rstart,
     *      air0,amat,l3,vbar,beta,a,b,c,j0,rhs,bc,imax,air1,
     *      switch,cstart,maxv,show)
c     formulate and solve boundary condition matrix
      implicit real*8(a-h,o-z)
      parameter (s1=21.0d0, s2=770.0d0, s3=50666.0d0, pi=3.1415926536d0)
      real*8 ak2(mch),coupl(mch,mch,kmax),
     *    fcc(mch,mch,2),gcc(mch,mch,2),alpha,rstart(nsteps+1),
     *    air0(nch,4),amat(l3),vbar(nch),beta(nch),a(nch,nch),
     *    b(nch,nch),c(nch,nch),rhs(imax),bc(12),air1(nch,2)
      integer j0(imax),cstart,show
      logical correc,deriv,airy
c
      call jarray(j0,nch,imax)
      airy = .true.
c
c************************************** large loop over each radial step
      do 100 ns=1,nsteps
         r0 = rstart(ns)
         r1 = rstart(ns+1)
         rh = r1-r0
         rb = r0 + 0.5d0*rh
         if (r0.ge.switch) airy = .false.
c------------------------------------- lsq part 1
         rls = rh/20.0d0
         sum = s1
         sumxx = s2*rls**2
         sum4x = s3*rls**4
         ted= 1.0d0/(sum*sum4x- sumxx**2)
         d11 = sum4x*ted
         d13 = -sumxx*ted
         d22 = 1.0d0/sumxx
         d33 = sum*ted
c      ted=1/determinant
c---------------------------------- dvbar---------
         if (airy) then
c           dvbar is gradient at mid-point
            dvbar = 0.0d0
            do 10 nc=1,nch
               do 10 nk=1,kmax
 10              dvbar= dvbar- dble(nk)*coupl(nc,nc,nk)
     *                 /(rb**(nk+1))
            dvbar = dvbar/nch
            if(dvbar.lt.0.0d0) then
               alpha = -1.0d0*(abs(dvbar)**(1.0d0/3.0d0))
            else
               alpha = dvbar**(1.0d0/3.0d0)
            endif
         endif
c
c---------------------------------------vbar------
         do 30 nc=1,nch
c           vbar is average potential
            vbar(nc) = coupl(nc,nc,1)*log(r1/r0)
            do 15 nk= 2,kmax
               vbar(nc) = vbar(nc) + coupl(nc,nc,nk)/dble(nk-1)*
     *                   (r0**(1-nk) - r1**(1-nk))
 15         continue
            vbar(nc) =  vbar(nc)/rh
            if (airy) then
              beta(nc) = (vbar(nc)-ak2(nc))/dvbar - rb
            else
             temp = ak2(nc) - vbar(nc)
              if (temp.lt.0.0) then
                 write(6,*)'e < v at switching radius with r>switch'
                 write(6,*)'use airy functions (ie increase switch)'
                 stop
              endif
              beta(nc) = sqrt(temp)
              alpha =1.0d0
            endif
c
c--------------------------------------  summing v-vo
            do 25 nc2= 1,nch
               sumf = 0.0d0
               sumfx = 0.0d0
               sumfxx = 0.0d0
c      super-fast 21 point least-squares fit to u-uo
               do 22 nl=-10,10
                  vt = 0.0d0
                  rd = dble(nl)*rls
                  do 23 nk=1,kmax
23                  vt = vt + coupl(nc,nc2,nk)/((rb+rd)**nk)
                  if (nc.eq.nc2) then
                    vt = vt -vbar(nc)
                    if (airy) vt = vt -rd*dvbar
                  endif
                  sumf = sumf+ vt
                  sumfx = sumfx+ vt*rd
                  sumfxx = sumfxx + vt*rd**2
 22            continue
c-------------------------------------------- lsq part 2 -------
               at = d11*sumf + d13*sumfxx
               bt = d22*sumfx
               ct = d13*sumf + d33*sumfxx
               if (airy) then
                  a(nc,nc2) = at - bt*rb +ct*rb*rb
                  b(nc,nc2) = bt - (ct+ct)*rb
                  c(nc,nc2) = ct
               else
                  a(nc,nc2) = at
                  b(nc,nc2) = bt
                  c(nc,nc2) = ct
               endif
c
c---------------------------------------------------------------
 25         continue
 30      continue
c
         ia=(ns-2)*2*nch
         ja= ia+2*nch
         jb= ja + nch
         id= ia + nch
         ic= id + nch
         idc= ic + nch
c-------------------------------------calculate and store airy's
         do 40 nc=1,nch
c
            if (airy) then
               dx0 = alpha*(r0+beta(nc))
               dx1 = alpha*(r1+beta(nc))
               call dairy(dx0,ai0,ad0,bi0,bd0)
               call dairy(dx1,ai1,ad1,bi1,bd1)
            else
               n2pi = int(r0/(2.d0*pi))
               r0 = r0 - n2pi*2.0d0*pi
               r1 = r1 - n2pi*2.0d0*pi
               rb = (r1+r0)/2.d0
               dx0 = beta(nc)*r0
               dx1 = beta(nc)*r1
               ai0 = cos(dx0)
               ai1 = cos(dx1)
               bi0 = sin(dx0)
               bi1 = sin(dx1)
               ad0 = -beta(nc)*bi0
               ad1 = -beta(nc)*bi1
               bd0 =  beta(nc)*ai0
               bd1 =  beta(nc)*ai1
            endif
c
            if (ia.lt.0) then
c              values at rmin not needed in the matrix
               air0(nc,1) = ai0
               air0(nc,2) = bi0
               air0(nc,3) = ad0
               air0(nc,4) = bd0
               ja = 0
               jb = nch
            else
               amat(j0(ia+nc)+ja+nc) =-ai0
               amat(j0(ia+nc)+jb+nc) =-bi0
               amat(j0(id+nc)+ja+nc) =-ad0
               amat(j0(id+nc)+jb+nc) =-bd0
            endif
            amat(j0(ic+nc)+ja+nc) = ai1
            amat(j0(ic+nc)+jb+nc) = bi1
            amat(j0(idc+nc)+ja+nc) = ad1
            amat(j0(idc+nc)+jb+nc) = bd1
c
            if (.not.deriv) then
            if (ns.eq.nstepr2) then
               air1(nc,1) = ai0
               air1(nc,2) = bi0
            endif
            endif
 40      continue
c
         if (airy) pal = pi/alpha
c----------------------------------  correction loop
         do 50 nc=1,nch
c
            if (ia.lt.0) then
               ai0 = air0(nc,1)
               bi0 = air0(nc,2)
               ad0 = air0(nc,3)
               bd0 = air0(nc,4)
            else
               ai0 =-amat(j0(ia+nc)+ja+nc)
               bi0 =-amat(j0(ia+nc)+jb+nc)
               ad0 =-amat(j0(id+nc)+ja+nc)
               bd0 =-amat(j0(id+nc)+jb+nc)
            endif
            ai1 = amat(j0(ic+nc)+ja+nc)
            bi1 = amat(j0(ic+nc)+jb+nc)
            ad1 = amat(j0(idc+nc)+ja+nc)
            bd1 = amat(j0(idc+nc)+jb+nc)
c
            bet1 = beta(nc)
           if (correc) then
c---------------------------------------- do off diagonal's first
              do 45 nc2=nc+1,nch
                 bet2 = beta(nc2)
                 if (ia.lt.0) then
                    ai2 = air0(nc2,1)
                    bi2 = air0(nc2,2)
                    ad2 = air0(nc2,3)
                    bd2 = air0(nc2,4)
                 else
                    ai2 =-amat(j0(ia+nc2)+ja+nc2)
                    bi2 =-amat(j0(ia+nc2)+jb+nc2)
                    ad2 =-amat(j0(id+nc2)+ja+nc2)
                    bd2 =-amat(j0(id+nc2)+jb+nc2)
                 endif
                 ai3 = amat(j0(ic+nc2)+ja+nc2)
                 bi3 = amat(j0(ic+nc2)+jb+nc2)
                 ad3 = amat(j0(idc+nc2)+ja+nc2)
                 bd3 = amat(j0(idc+nc2)+jb+nc2)
c
             if (airy) then
               call ofdag(r0,r1,alpha,bet1,bet2,bc,ai0,bi0,ad0,bd0,
     *           ai1,bi1,ad1,bd1,ai2,bi2,ad2,bd2,ai3,bi3,ad3,bd3)
             else
               call scoff(r0,r1,bet1,bet2,bc,ai0,bi0,ai1,bi1,
     *           ai2,bi2,ai3,bi3)
             endif
c
                 w = a(nc,nc2)
                 wr= b(nc,nc2)
                 wrr = c(nc,nc2)
                 if (.not.airy) then
                    w = w - wr*rb +wrr*rb*rb
                    wr = wr - 2.d0*wrr*rb
                    pal = 1.d0/bet1
                 endif
      wp371  =    pal*(w*bc(3)+wr*bc(7)+wrr*bc(11))
      wp159  =    pal*(w*bc(1)+wr*bc(5)+wrr*bc(9))
      wp261  =    pal*(w*bc(2)+wr*bc(6)+wrr*bc(10))
      wp481  =    pal*(w*bc(4)+wr*bc(8)+wrr*bc(12))
c
c-------------------------------------correct (nc,nc2)
      amat(j0(ic+nc)+ja+nc2)= -ai1*wp371+ bi1*wp159
      amat(j0(idc+nc)+ja+nc2)=(-ad1*wp371+ bd1*wp159)*alpha
      amat(j0(ic+nc)+jb+nc2)=  bi1*wp261- ai1*wp481
      amat(j0(idc+nc)+jb+nc2)= (bd1*wp261- ad1*wp481)*alpha
c
                 w = a(nc2,nc)
                 wr= b(nc2,nc)
                 wrr = c(nc2,nc)
                 if (.not.airy) then
                    w = w - wr*rb +wrr*rb*rb
                    wr = wr - 2.d0*wrr*rb
                    pal = 1.d0/bet2
                 endif
c--------------------------------------correct (nc2,nc)
      wp371  =    pal*(w*bc(3)+wr*bc(7)+wrr*bc(11))
      wp159  =    pal*(w*bc(1)+wr*bc(5)+wrr*bc(9))
      wp261  =    pal*(w*bc(2)+wr*bc(6)+wrr*bc(10))
      wp481  =    pal*(w*bc(4)+wr*bc(8)+wrr*bc(12))
c
c-------------------------------------correct (nc,nc2)
      amat(j0(ic+nc2)+ja+nc)= -ai3*wp261+ bi3*wp159
      amat(j0(idc+nc2)+ja+nc)=(-ad3*wp261+ bd3*wp159)*alpha
      amat(j0(ic+nc2)+jb+nc)=  bi3*wp371- ai3*wp481
      amat(j0(idc+nc2)+jb+nc)=(bd3*wp371- ad3*wp481)*alpha
c
 45           continue
c-------------------------------- diagonal corrections
           endif
            if (.not.correc) then
                do 46 ji=1,12
 46                bc(ji)=0.0d0
                pal = 0.0d0
            else
               if (airy) then
                 call cdiag(r0,r1,alpha,bet1,bc,ai0,bi0,ad0,bd0,
     *            ai1,bi1,ad1,bd1)
               else
                 call scorr(r0,r1,bet1,bc,ai0,bi0,ai1,bi1)
               endif
c
            endif
            w = a(nc,nc)
            wr= b(nc,nc)
            wrr = c(nc,nc)
            if (.not.airy) then
               w = w - wr*rb +wrr*rb*rb
               wr = wr - 2.d0*wrr*rb
               pal = 1.d0/bet1
            endif
c
            wp258 = pal*(w*bc(2)+wr*bc(5)+wrr*bc(8))
            wp147 = pal*(w*bc(1)+wr*bc(4)+wrr*bc(7))
            wp369 = pal*(w*bc(3)+wr*bc(6)+wrr*bc(9))
c
      amat(j0(ic+nc)+ja+nc)= ai1*(1.0d0-wp258)+bi1*wp147
      amat(j0(idc+nc)+ja+nc)=(ad1*(1.0d0-wp258) +bd1*wp147)*alpha
      amat(j0(ic+nc)+jb+nc)= bi1*(1.0d0+wp258) -ai1*wp369
      amat(j0(idc+nc)+jb+nc)=(bd1*(1.0d0+wp258) -ad1*wp369)*alpha
c
c----------------------------------- mult. deriv.s by alpha
          if (airy) then
            if (ia.lt.0) then
               air0(nc,3) = air0(nc,3)*alpha
               air0(nc,4) = air0(nc,4)*alpha
            else
               amat(j0(id+nc)+ja+nc) = amat(j0(id+nc)+ja+nc)*alpha
               amat(j0(id+nc)+jb+nc) = amat(j0(id+nc)+jb+nc)*alpha
            endif
          endif
c
 50      continue
c--------------------------end of correction loop
 100  continue
c********************************************* end of nsteps loop
c
      call ludcst(amat,maxv,imax,imax,j0,nch)
c
c
c                     jst is awkward off-set for partition.ne.1
      jst = cstart-1
      ist = imax-2*nch
c----------------------------------- backsubstitute
      do 140 nc=1,nch
         do 150 i=1,imax
 150        rhs(i) = 0.0d0
          if(show.ge.6) write(6,49) '0',nc,(air0(nc,i),i=1,4)
          if(show.ge.6) write(6,49) '1',nc,(air1(nc,i),i=1,4)
 49       format('air',a1,'(1-4) @',i3,':',1p,4g12.4)
         rhs(ist+nc) = fcc(nc+jst,nc+jst,1)
         rhs(ist+nc+nch) = fcc(nc+jst,nc+jst,2)
           if(show.ge.5) write(6,151) nc,(rhs(i),i=1,imax)
 151       format('rhs @',i3,':',8g12.4/(9x,10g12.4))
         call lubkst(amat,maxv,imax,imax,j0,nch,rhs)
           if(show.ge.5) write(6,152) nc,(rhs(i),i=1,imax)
 152       format('sol @',i3,':',8g12.4/(9x,10g12.4))
         do 160 nc2=1,nch
           fcc(nc2+jst,nc+jst,1)= rhs(nc2)*air0(nc2,1)+
     *                            rhs(nch+nc2)*air0(nc2,2)
           if (deriv) then
                fcc(nc2+jst,nc+jst,2)= rhs(nc2)*air0(nc2,3)
     *                        +rhs(nch+nc2)*air0(nc2,4)
           else
                ia = 2*(nstepr2-1)*nch
c               if (ia.gt.0)then
                  fcc(nc2+jst,nc+jst,2)= rhs(ia+nc2)*air1(nc2,1)
     *                          +rhs(ia+nch+nc2)*air1(nc2,2)
c               else
c                 fcc(nc2+jst,nc+jst,2)=fcc(nc2+jst,nc+jst,1)
c               endif
           endif
        if(show.ge.2.and.abs(fcc(nc2,nc,1)).gt.1e-9.or.show.ge.3)
     x   write(6,39) 'f',nc2,nc,ia,(fcc(nc2,nc,i),i=1,2)
     x            ,rhs(ia+nc2),air1(nc2,1) ,rhs(ia+nch+nc2),air1(nc2,2)
 39       format('matrix ',a1,' to',i3,' from ch.',i3,i3,':',1p,6g12.4)
 160     continue
c---------------------------------done the f's now do the g's
         do 170 i=1,imax
 170        rhs(i) = 0.0d0
         rhs(ist+nc) = gcc(nc+jst,nc+jst,1)
         rhs(ist+nc+nch) = gcc(nc+jst,nc+jst,2)
         call lubkst(amat,maxv,imax,imax,j0,nch,rhs)
         do 180 nc2=1,nch
            gcc(nc2+jst,nc+jst,1)= rhs(nc2)*air0(nc2,1)+
     *                             rhs(nch+nc2)*air0(nc2,2)
            if (deriv) then
                gcc(nc2+jst,nc+jst,2)= rhs(nc2)*air0(nc2,3)
     *                                +rhs(nch+nc2)*air0(nc2,4)
            else
c               if (ia.gt.0) then
                   gcc(nc2+jst,nc+jst,2)= rhs(ia+nc2)*air1(nc2,1)
     *                           +rhs(ia+nch+nc2)*air1(nc2,2)
c               else
c                  gcc(nc2+jst,nc+jst,2)= gcc(nc2+jst,nc+jst,1)
c               endif
            endif
        if(show.ge.2.and.abs(gcc(nc2,nc,1)).gt.1e-9.or.show.ge.3)
     x   write(6,39) 'g',nc2,nc,ia,(gcc(nc2,nc,i),i=1,2)
 180     continue
c
 140  continue
c--------------------------------end of story!
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine coul(mch,nch,ak2,fcc,gcc,rout,
     *      f,fp,g,gp,eta,lval,lmax1,coupl,kmax)
c    get uncoupled coulomb wfns - better to do in main program.
      implicit real*8(a-h,o-z)
      real*8 eta(mch),ak2(mch),f(lmax1),rout,coupl(mch,mch,kmax),
     *    fcc(mch,mch,2),gcc(mch,mch,2),g(lmax1),fp(lmax1),gp(lmax1)
      integer lval(mch)
c
c
      xlm =  0.0d0
      xlx =  dble(lmax1-1)
      do 100 nc =1,nch
        ak = sqrt(ak2(nc))
	eta(nc)=coupl(nc,nc,1)/(2.d0*ak)
	al1= 0.5d0*(sqrt(4.d0*coupl(nc,nc,2)+1)-1.d0)+0.001d0
	lval(nc) = int(al1)
        rho = rout*ak
        kfn=0
        mode=1
        ifail=0
        call coulfg(rho,eta(nc),xlm,xlx,f,g,fp,gp,mode,kfn,ifail)
        i = lval(nc)+1
        fcc(nc,nc,1) =f(i)
        fcc(nc,nc,2) =fp(i)*ak
        gcc(nc,nc,1) =g(i)
        gcc(nc,nc,2) =gp(i)*ak
100     continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cdiag(r0,r1,alp,bet,bc,ai0,bi0,ad0,bd0,ai1,bi1,
     *      ad1,bd1)
c     diagonal correction integralsfor airy fns
      implicit real*8(a-h,o-z)
      parameter (pi=3.141592653589793d0)
c
      real*8 alp,bet,bet2,bet3,pha2,pha3,ai0,ad0,bi0,bd0,
     +       ai1,ad1,bi1,bd1,c11,c12,c13,c01,c02,c03,cb,cc,
     +       pha,bc,rb0,rb1
      real*8 r0,r02,r03,r1,r12,r13
      dimension bc(12)
c
c     correction terms
      bet2 = bet*bet
      bet3 = bet2*bet
      pha = 1.0d0/alp
      pha2 = pha*pha
      pha3 = pha2*pha
      rb0 = r0+bet
      rb1 = r1+bet
      r02 = r0*r0
      r03 = r0*r02
      r12 = r1*r1
      r13 = r1*r12
      ai12 = ai1*ai1
      ai02 = ai0*ai0
      ad12 = ad1*ad1
      ad02 = ad0*ad0
      aid1 = ai1*ad1
      aid0 = ai0*ad0
      aibi1= ai1*bi1
      aibi0= ai0*bi0
      aibd1= ai1*bd1
      aibd0= ai0*bd0
      adbi1= bi1*ad1
      adbi0= bi0*ad0
      adbd1 = ad1*bd1
      adbd0 = ad0*bd0
      bi12 = bi1*bi1
      bi02 = bi0*bi0
      bid1 = bi1*bd1
      bid0 = bi0*bd0
      bd12 = bd1*bd1
      bd02 = bd0*bd0
c
      c11=(r1+bet)
      c13=- pha
      c01=(r0+bet)
      c03= c13
c-----------------  int(aa), int(ab), int(bb)
      bc(1) = c11*ai12 + c13*ad12 - c01*ai02 - c03*ad02
      bc(2) = c11*aibi1 + c13*adbd1 - c01*aibi0 - c03*adbd0
      bc(3) = c11*bi12 + c13*bd12 - c01*bi02 - c03*bd02
c
      cb = 1.0d0/3.0d0
c     c11= cb*(rb1*rb1-3.0d0*rb1*bet)
      c11= cb*(r12-bet*r1-bet2-bet2)
      c12= cb*pha2/2.0d0
      c13= cb*(bet+bet-r1)*pha
      c01= cb*(r02-bet*r0-bet2-bet2)
      c02= c12
      c03= cb*(bet+bet-r0)*pha
c
c------------------ int(raa), int(rab), int(rbb)
      bc(4) = (c13*ad12-c03*ad02) + (c11*ai12-c01*ai02)
     *      + 2.0d0 *(c12*aid1-c02*aid0)
      bc(5) = (c13*adbd1-c03*adbd0) + (c11*aibi1-c01*aibi0)
     *      + c12*(adbi1+aibd1) - c02*(adbi0+aibd0)
      bc(6) = (c13*bd12-c03*bd02) + (c11*bi12-c01*bi02)
     *      + 2.0d0 *(c12*bid1-c02*bid0)
c
      cc = 1.0d0/15.0d0
      b3p3 = 8.0d0*bet3 - 3.0d0*pha3
      c11= cc*(3.0d0*r13 -r12*bet +4.0d0*r1*bet2 + b3p3)
      c12= cc*(3.0d0*r1 -2.0d0*bet)*pha2
      c13= cc*(4.0d0*r1*bet -3.0d0*r12 -8.0d0*bet2)*pha
      c01= cc*(3.0d0*r03 -r02*bet +4.0d0*r0*bet2 + b3p3)
      c02= cc*(3.0d0*r0 -2.0d0*bet)*pha2
      c03= cc*(4.0d0*r0*bet -3.0d0*r02 -8.0d0*bet2)*pha
c------------------ int(rraa), int(rrab), int(rrbb)
c
      bc(7)= (c13*ad12-c03*ad02) +(c11*ai12-c01*ai02)
     *     + 2.0d0* (c12*aid1-c02*aid0)
      bc(8)= (c13*adbd1-c03*adbd0) +(c11*aibi1-c01*aibi0)
     *     + c12*(aibd1+adbi1) - c02*(aibd0+adbi0)
      bc(9)= (c13*bd12-c03*bd02) +(c11*bi12-c01*bi02)
     *     + 2.0d0* (c12*bid1-c02*bid0)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ofdag(r0,r1,alp,bet1,bet2,bc,ai0,bi0,ad0,bd0,
     *             ai1,bi1,ad1,bd1,ai2,bi2,ad2,bd2,ai3,bi3,ad3,bd3)
c     off diagonal couplingsfor airy fns
      implicit real*8(a-h,o-z)
      parameter (pi=3.141592653589793d0)
c
      real*8 alp,bet2,pha2,pha3,ai0,ad0,bi0,bd0,
     +       ai1,ad1,bi1,bd1,c11,c12,c13,c01,c02,c03,cb,cc,
     +       pha,bc
      real*8 r0,r1
      dimension bc(12)
c
c     correction terms
      bpb = bet1+bet2
      bmb  = bet1-bet2
c
c     if degenerate
c
      if (abs(bmb/bet1).lt.1.0d-3) then
 10       call cdiag(r0,r1,alp,bet1,bc,ai0,bi0,ad0,bd0,
     *       ai1,bi1,ad1,bd1)
          bc(12)=bc(9)
          bc(11)=bc(8)
          bc(10)=bc(8)
          bc(9)=bc(7)
          bc(8)=bc(6)
          bc(7)=bc(5)
          bc(6)=bc(5)
          bc(5)=bc(4)
          bc(4)=bc(3)
          bc(3)=bc(2)
c         write(6,*) bc(9)
          return
      endif
c
      bmb2 = bmb*bmb
      dbmb = 1.0d0/bmb
      dbmb2 = dbmb*dbmb
      pha = 1.0d0/alp
      pha2 = pha*pha
      pha3 = pha2*pha
      p2db = pha2*dbmb
      p3db2 = pha3*dbmb2
      adi13 = ad1*ai3 - ai1*ad3
      adi02 = ad0*ai2 - ai0*ad2
      ab13  = ad1*bi3 - ai1*bd3
      ab02  = ad0*bi2 - ai0*bd2
      ba13  = bd1*ai3 - bi1*ad3
      ba02  = bd0*ai2 - bi0*ad2
      bdi13 = bd1*bi3 - bi1*bd3
      bdi02 = bd0*bi2 - bi0*bd2
c------------------------- int(aa), int(ab), int(ba), int(bb)
      c12= p2db
      bc(1) = c12*(adi13-adi02)
      bc(2) = c12*(ab13-ab02)
      bc(3) = c12*(ba13-ba02)
      bc(4) = c12*(bdi13-bdi02)
c
      p5d3 = p3db2*p2db
      c11 = -(bpb+r1+r1)*p3db2
      c01 = -(bpb+r0+r0)*p3db2
      c12 = r1*p2db + p5d3 +p5d3
      c02 = r0*p2db + p5d3 +p5d3
      c13 = 2.0d0*pha *p3db2
      c03 = c13
c------------------------- int(raa), int(rab), int(rba), int(rbb)
      bc(5) = c11*ai1*ai3 + c12*adi13 + c13*ad1*ad3 -
     *        c01*ai0*ai2 - c02*adi02 - c03*ad0*ad2
      bc(6) = c11*ai1*bi3 + c12*ab13  + c13*ad1*bd3 -
     *        c01*ai0*bi2 - c02*ab02  - c03*ad0*bd2
      bc(7) = c11*bi1*ai3 + c12*ba13  + c13*bd1*ad3 -
     *        c01*bi0*ai2 - c02*ba02  - c03*bd0*ad2
      bc(8) = c11*bi1*bi3 + c12*bdi13 + c13*bd1*bd3 -
     *        c01*bi0*bi2 - c02*bdi02 - c03*bd0*bd2
c
      alp3 = alp**3
      a3b2 = alp3*bmb2
      p6d4 = p3db2*p3db2
      db5 = dbmb2*dbmb2*dbmb *pha3*pha3*pha2
      cc = 24.0d0 + 2.0d0*alp3*bpb*bmb2
      c11 = -p6d4 * (12.0d0*bpb + cc*r1 +
     *       4.0d0*a3b2*r1*r1)
      c01 = -p6d4 * (12.0d0*bpb + cc*r0 +
     *       4.0d0*a3b2*r0*r0)
      c12 = db5 * (cc+
     *      12.0d0*a3b2*r1 + a3b2*a3b2*r1*r1)
      c02 = db5 * (cc+
     *      12.0d0*a3b2*r0 + a3b2*a3b2*r0*r0)
      c13 = 4.0d0*pha*p6d4*(6.0d0+a3b2*r1)
      c03 = 4.0d0*pha*p6d4*(6.0d0+a3b2*r0)
      cb = 4.0d0*p5d3
c
c     test for loss of precision: if so do diagonal coupling routine
c
      if (abs((c13-c03)/c13).lt.1.d-6) then
           write(6,*)'rounding error so using diagonal integral'
           goto 10
      endif
c
c------------------------- int(rraa), int(rrab), int(rrba), int(rrbb)
      bc(9) = c11*ai1*ai3 + c12*adi13 + c13*ad1*ad3 -
     *        c01*ai0*ai2 - c02*adi02 - c03*ad0*ad2 +
     *        cb*(bet2*(ad1*ai3-ad0*ai2) - bet1*(ai1*ad3-ai0*ad2))
      bc(10)= c11*ai1*bi3 + c12*ab13  + c13*ad1*bd3 -
     *        c01*ai0*bi2 - c02*ab02  - c03*ad0*bd2 +
     *        cb*(bet2*(ad1*bi3-ad0*bi2) - bet1*(ai1*bd3-ai0*bd2))
      bc(11)= c11*bi1*ai3 + c12*ba13  + c13*bd1*ad3 -
     *        c01*bi0*ai2 - c02*ba02  - c03*bd0*ad2 +
     *        cb*(bet2*(bd1*ai3-bd0*ai2) - bet1*(bi1*ad3-bi0*ad2))
      bc(12)= c11*bi1*bi3 + c12*bdi13 + c13*bd1*bd3 -
     *        c01*bi0*bi2 - c02*bdi02 - c03*bd0*bd2 +
     *        cb*(bet2*(bd1*bi3-bd0*bi2) - bet1*(bi1*bd3-bi0*bd2))
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine jarray(j0,nch,imax)
       integer j0(imax)
c      new array technique for staircase matrix.
c
       do 10 i=1,imax
          k= (i-1)/(2*nch)
          j0(i) = (i-1)*4*nch - 2*k*nch
 10       continue
c
       i=imax-2*nch+2
       jsub =2*nch
c
 20    j0(i) = j0(i) -jsub
       i=i+1
       jsub = jsub + 2*nch
       if (i.le.imax) goto 20
c
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine ludcst(a,np1,np2,n,j0,nch)
c
c given an n * n matrix a, with physical dimension np, this routine
c replaces it with the lu decomposition of a rowwise permutation of
c itself. staircase modifications as 1-d array.
c no pivoting!
      implicit real*8(a-h,o-z),integer(i-n)
      parameter (tiny=1.0d-50)
      dimension a(np1),j0(np2)
      do 19 j=1,n
         jmd = ((j-1)/(2*nch)-1)*2*nch
         do 14 i=max(1,1+jmd),j-1
         sum = a(j0(i)+j)
           imd = ((i-1)/(2*nch))*2*nch
           do 13 k=max(1,imd+1,jmd+1),
     *             min(i-1,imd+4*nch,jmd+4*nch)
13         sum = sum - a(j0(i)+k)*a(j0(k)+j)
14       a(j0(i)+j) = sum
c
         do 16 i=j,min(n,jmd+4*nch)
            sum = a(j0(i)+j)
            imd = ((i-1)/(2*nch))*2*nch
            do 15 k=max(1,imd+1,jmd+1),
     *              min(j-1,imd+4*nch,jmd+4*nch)
15            sum = sum - a(j0(i)+k)*a(j0(k)+j)
            a(j0(i)+j) = sum
c
16          continue
c
         if(abs(a(j0(j)+j)).lt.tiny) a(j0(j)+j) = tiny
         if(j.ne.n) then
            dum = 1.0d0 /a(j0(j)+j)
            do 18 i=j+1,min(n,jmd+4*nch)
18             a(j0(i)+j) = a(j0(i)+j)*dum
         endif
19     continue
      return
      end
      subroutine lubkst(a,np1,np2,n,j0,nch,b)
c
c solves the set of n linear equations a.x = b.
c here, a is input, not as the matrix a, but rather its lu decomposition
c
      implicit real*8(a-h,o-z),integer(i-n)
      dimension a(np1),b(np2),j0(np2)
      ii = 0
      do 12 i=1,n
        sum = b(i)
        if(ii.ne.0) then
          imd = ((i-1)/(2*nch))*2*nch
          do 11 j=imd+1,i-1
11           sum = sum - a(j0(i)+j)*b(j)
          else if(abs(sum).gt.0.0) then
             ii = i
          endif
12        b(i) = sum
      do 14 i=n,1,-1
        sum = b(i)
        if(i.lt.n) then
           imd = ((i-1)/(2*nch))*2*nch
           do 13 j = i+1,min(n,imd+4*nch)
13            sum = sum - a(j0(i)+j)*b(j)
           endif
14      b(i) = sum/a(j0(i)+i)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scorr(r0,r1,bet1,bc,ai0,bi0,ai1,bi1)
c     diagonal correction integrals for sines + cosines
      implicit real*8(a-h,o-z)
      real*8 r0,r1,bet1,bc(12),ai0,ai1,bi0,bi1
c
      sn20 = 2.d0*ai0*bi0
      sn21 = 2.d0*ai1*bi1
      cs20 = ai0*ai0 - bi0*bi0
      cs21 = ai1*ai1 - bi1*bi1
      r02 = r0 *r0
      r12 = r1 *r1
      bte = 1.d0/bet1
      sms = 0.5d0*bte*(sn21-sn20)
      cmc = 0.5d0*bte*(cs21-cs20)
      bt4 = 0.25d0*bte
c
      bc(1) = 0.5d0*(r1-r0 + bte*(ai1*bi1 - ai0*bi0))
      bc(2) = 0.5d0*bte * (ai0*ai0 - ai1*ai1)
      bc(3) = 0.5d0*(r1-r0 - bte*(ai1*bi1 - ai0*bi0))
c
      cr3 = 0.25d0*(r12-r02)
      cr4 = bt4*(r1*sn21-r0*sn20 +cmc)
      bc(4) = cr3 + cr4
      bc(6) = cr3 - cr4
      bc(5) = bt4*(r0*cs20-r1*cs21 +sms)
c
      cr3 = (r1*r12-r0*r02)/6.0d0
      cr4 =  bt4*(r12*sn21-r02*sn20 + bte*
     *    (r1*cs21 - r0*cs20  - sms))
      bc(7) = cr3 + cr4
      bc(9) = cr3 - cr4
      bc(8)= bt4*(r02*cs20-r12*cs21 +bte*
     *    (r1*sn21 - r0*sn20 + cmc))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scoff(r0,r1,bet1,bet2,bc,ai0,bi0,ai1,bi1,
     *     ai2,bi2,ai3,bi3)
c        off diagonal sines+cosines correction integrals
      implicit real*8(a-h,o-z)
      real*8 bc(12)
c
      apb = bet1+bet2
      amb = bet1-bet2
      bma = bet2-bet1
c
c     if degenerate then scorr
c
      if (abs(amb/bet1).lt.1.0d-6) then
        call scorr(r0,r1,bet1,bc,ai0,bi0,ai1,bi1)
        bc(12)= bc(9)
        bc(11)= bc(8)
        bc(10)= bc(8)
        bc(9) = bc(7)
        bc(8) = bc(6)
        bc(7) = bc(5)
        bc(6) = bc(5)
        bc(5) = bc(4)
        bc(4) = bc(3)
        bc(3) = bc(2)
c        write(6,*)bc(9)
        return
      endif
c         a=cos   b=sin  1,3 = r1   0,2 = r0
      sapb1 = (ai1*bi3 + bi1*ai3)/apb
      samb1 = (bi1*ai3 - ai1*bi3)/amb
      capb1 = (ai1*ai3 - bi1*bi3)/apb
      camb1 = (ai1*ai3 + bi1*bi3)/amb
      sbma1 = samb1
      cbma1 = -camb1
c
      sapb0 = (ai0*bi2 + bi0*ai2)/apb
      samb0 = (bi0*ai2 - ai0*bi2)/amb
      capb0 = (ai0*ai2 - bi0*bi2)/apb
      camb0 = (ai0*ai2 + bi0*bi2)/amb
      sbma0 = samb0
      cbma0 = -camb0
c
      sp1 = 0.5d0*(sapb1+samb1)
      sp0 = 0.5d0*(sapb0+samb0)
      sm1 = 0.5d0*(samb1-sapb1)
      sm0 = 0.5d0*(samb0-sapb0)
      ca1 = 0.5d0*(capb1+camb1)
      ca0 = 0.5d0*(capb0+camb0)
      cb1 = 0.5d0*(cbma1+capb1)
      cb0 = 0.5d0*(cbma0+capb0)
 
      bc(1) = sp1 - sp0
      bc(2) = cb0 - cb1
      bc(3) = ca0 - ca1
      bc(4) = sm1 - sm0
c
      cp = (capb1 - capb0)/apb
      ca = (camb1 - camb0)/amb
      cb = (cbma1 - cbma0)/bma
      sp = (sapb1 - sapb0)/apb
      sa = (samb1 - samb0)/amb
      sb = (sbma1 - sbma0)/bma
c
      bc(5) = r1*sp1 - r0*sp0 + 0.5d0*(ca+cp)
      bc(6) = r0*cb0 - r1*cb1 + 0.5d0*(sb+sp)
      bc(7) = r0*ca0 - r1*ca1 + 0.5d0*(sa+sp)
      bc(8) = r1*sm1 - r0*sm0 + 0.5d0*(ca-cp)
c
      cp = cp/apb
      sp = sp/apb
      sa = sa/amb
      rcam = (r1*camb1-r0*camb0)/amb
      rcap = (r1*capb1-r0*capb0)/apb
      rsap = (r1*sapb1-r0*sapb0)/apb
      r02 = r0*r0
      r12 = r1*r1
c
      bc(9) = r12*sp1 - r02*sp0 +rcam + rcap -sa -sp
      bc(10)= r02*cb0 - r12*cb1 +rsap +cb/bma + cp +
     *   (r1*sbma1-r0*sbma0)/bma
      bc(11)= r02*ca0 - r12*ca1 +rsap +ca/amb + cp +
     *   (r1*samb1-r0*samb0)/amb
      bc(12) = r12*sm1 - r02*sm0 +rcam - rcap -sa +sp
c
      return
      end
c
      subroutine dairy(dx,ai,aip,bi,bip)
c  for double precision arguments, this routine calculates the airy
c     function ai(x) and its derivative aip(x).  it also finds
c     the other real linearly independent solution bi(x) and
c     its derivative bip(x).
c     the definitions and normalizations are as in nbs handbook
c     of mathematical functions,p.446
c     the methods used are power series expansion for small x
c     and gaussian integration for large x
      dimension x(16),w(16),xsq(16)
      double precision dx,ai,aip,bi,bip
      double precision    xs ,xcube,aisum,aipsum
      double precision df,dfp,dg,dgp
      double precision fjm2,fjm1,fj,fjp1,fjp2,factor
      double precision c1,c2,root3
      double precision dzeta,darg,drootx
      double precision root4x,s,co,ratio,efac,zetasq
      double precision sumr,sumi,sumrp,sumip,termr,termi
      double precision dzero,da,db,den,one
      double precision x,w,xsq
      double precision rsq,temp,   rtpi,rtpi2
      double precision terma,termb
      logical needbi
      data dzero,one /0.0d0,1.0d0/
      data root3/1.732050807568877d0/
      data c1,c2 /.355028053887817d0, .258819403792807d0/
      data rtpi /.2820947917738781d0/
      data rtpi2/.5641895835477562d0/
      data rsq,den,termi/3*0.0d0/
c *** inserted arb 6.80
c  positions and weights for 10-term sum for airy functions
       data w( 1) /  3.1542515762964787d-14/
       data w( 2) /  6.6394210819584921d-11/
       data w( 3) /  1.7583889061345669d-08/
       data w( 4) /  1.3712392370435815d-06/
       data w( 5) /  4.4350966639284350d-05/
       data w( 6) /  7.1555010917718255d-04/
       data w( 7) /  6.4889566103335381d-03/
       data w( 8) /  3.6440415875773282d-02/
       data w( 9) /  1.4399792418590999d-01/
       data w(10) /  8.1231141336261486d-01/
       data x( 1) /  1.4083081072180964d+01/
       data x( 2) /  1.0214885479197331d+01/
       data x( 3) /  7.4416018450450930d+00/
       data x( 4) /  5.3070943061781927d+00/
       data x( 5) /  3.6340135029132462d+00/
       data x( 6) /  2.3310652303052450d+00/
       data x( 7) /  1.3447970824609268d+00/
       data x( 8) /  6.4188858369567296d-01/
       data x( 9) /  2.0100345998121046d-01/
       data x(10) /  8.0594359172052833d-03/
      data xsq( 1) /0.19833317248562170d 03/
      data xsq( 2) /0.10434388535311650d 03/
      data xsq( 3) /0.55377438020178170d 02/
      data xsq( 4) /0.28165249974668990d 02/
      data xsq( 5) /0.13206054139355800d 02/
      data xsq( 6) /0.54338651079380440d 01/
      data xsq( 7) /0.18084791929954200d 01/
      data xsq( 8) /0.41202095387883690d 00/
      data xsq( 9) /0.40402390924418070d-01/
      data xsq(10) /0.64954507303538390d-04/
c  positions and weights for  4-term sum for airy functions
       data w(11) /  4.7763903057577263d-05/
       data w(12) /  4.9914306432910959d-03/
       data w(13) /  8.6169846993840312d-02/
       data w(14) /  9.0879095845981102d-01/
       data x(11) /  3.9198329554455091d+00/
       data x(12) /  1.6915619004823504d+00/
       data x(13) /  5.0275532467263018d-01/
       data x(14) /  1.9247060562015692d-02/
      data xsq(11) /0.15365090398596670d 02/
      data xsq(12) /0.28613816631634610d 01/
      data xsq(13) /0.25276291648668180d 00/
      data xsq(14) /0.37044934027789980d-03/
c  positions and weights for  2-term sum for airy functions
       data w(15) /  9.6807280595773604d-01/
       data w(16) /  3.1927194042263958d-02/
       data x(15) /  3.6800601866153044d-02/
       data x(16) /  1.0592469382112378d+00/
      data xsq(15) /0.13542842977111070d-02/
      data xsq(16) /0.11220040761098810d 01/
      if(dx.lt.-5.0d0) go to 100
      needbi=.false.
      if(dx.gt.3.7d0) go to 200
c     this route for smallx, using power series.
c     initialize
10    xs  = dx*dx
      xcube = xs *dx
      xs  = xs *0.5d0
      df = c1
      dfp = c1*xs
      dg = c2*dx
      dgp = c2
      aisum = df - dg
      aipsum = dfp - dgp
      bi = df + dg
      bip = dfp + dgp
      fjm2=-2.0d0
20    fjm2=fjm2+3.0d0
      fjm1=fjm2+one
      fj=fjm1+one
      fjp1=fj+one
      fjp2=fjp1+one
      ratio = xcube/fj
      df = df*ratio/fjm1
      dfp = dfp*ratio/fjp2
      dg = dg*ratio/fjp1
      dgp = dgp*ratio/fjm2
      bi = bi + (df+dg)
      bip = bip + (dfp+dgp)
      if(needbi) go to 80
      aisum = aisum + (df-dg)
      aipsum = aipsum + (dfp-dgp)
c     convergence test
80    if(dabs(df).gt.1.0d-16) go to 20
c     convergence. compute functions
99    bi = root3*bi
      bip = root3*bip
c  this returns if x is between 3.7 and 8.0, since in such cases more
c  accurate values of ai and aip have already been found by gaussian
c  integration
      if(needbi)return
      ai = aisum
      aip = aipsum
      return
c  gaussian integration for large negative x
100   drootx = dsqrt(-dx)
      root4x = dsqrt(drootx)
      dzeta = -.6666666666666667d0*dx*drootx
      darg = dzeta - .7853981633974483d0
      sumr = dzero
      sumi = dzero
      sumrp = dzero
      sumip = dzero
c  test to see how many terms are needed in gaussian integration
      if(dx.lt.(-200.d0)) go to 140
      if(dx.lt.(-15.d0)) go to 130
c  this case for dx between -5.0 and -15.0
      limlo=1
      limhi=10
      go to 149
c  this case for dx between -15.0 and -200.
130   limlo=11
      limhi=14
      go to 149
c  this case for dx.lt.-200.
140   limlo=15
      limhi=16
149   zetasq=dzeta**2
      do 150 k=limlo,limhi
      termr=w(k)/((zetasq+xsq(k))**2)
      sumr = sumr + termr
      termr=termr*x(k)
      sumi=sumi+termr
      termr=termr*x(k)
      sumrp=sumrp+termr
150   sumip=sumip+termr*x(k)
      sumr=(sumr*zetasq+sumrp)*zetasq
      temp=sumi*zetasq
      sumi=(temp+sumip)*dzeta
      sumrp=sumrp*dzeta
      sumip=sumip-temp
c  form airy functions
196   s = dsin(darg)
      co = dcos(darg)
      ratio = rtpi2/root4x
      ai = ratio*(co*sumr + s*sumi)
      bi = ratio*(co*sumi - s*sumr)
      sumrp=sumrp+sumrp
      ratio = -.25d0/dx
      factor = -rtpi2*root4x
      aip = ratio*ai - drootx*bi + factor*(co*sumrp+s*sumip)
      bip = ratio*bi + drootx*ai + factor*(co*sumip-s*sumrp)
      return
c   gaussian integration for large positive x
200   drootx = dsqrt(dx)
      dzeta = .6666666666666667d0*dx*drootx
      efac = dexp(-dzeta)
      root4x = dsqrt(drootx)
      ai = dzero
      bi = dzero
      aip = dzero
      bip = dzero
      if(dx.lt.8.0d0) needbi=.true.
c  test to see how many terms are needed in gaussian integration
      if(dx.gt.15.0d0) go to 230
c  this case for dx between 3.7 and 15.
      limlo=1
      limhi=10
      go to 249
c  this case for dx greater than 15.
230   limlo=11
      limhi=14
249   do 250 k=limlo,limhi
      da=dzeta+x(k)
      terma = w(k)/da
      ai = ai + terma
      aip=aip+terma*x(k)/da
      if(needbi) go to 250
      db=dzeta-x(k)
      termb = w(k)/db
      bi = bi + termb
      bip=bip+termb*x(k)/db
250   continue
c  form functions
      factor=rtpi*dzeta/root4x
      ratio = 0.25d0/dx
      ai=ai*efac*factor
      aip=-(drootx+ratio)*ai+rtpi*root4x*efac*aip
c  this is satisfied only for x between 3.7 and 8.0  in these cases
c  the bi and bip about to be computed are not sufficiently accurate.
c  thus return to power series for bi and bip.
      if(needbi) go to 10
      factor=factor+factor
      bi=bi*factor/efac
      bip=(drootx-ratio)*bi-rtpi2*root4x*bip/efac
      return
      end
c
      subroutine coulfg(xx,eta1,xlmin,xlmax, fc,gc,fcp,gcp,             
     *                  mode1,kfn,ifail)                                
c                                                                       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  revised coulomb wavefunction program using steed's method           c
c                                                                      c
c  a. r. barnett           manchester  march   1981                    c
c                                                                      c
c  original program 'rcwfn'      in    cpc  8 (1974) 377-395           c
c                 + 'rcwff'      in    cpc 11 (1976) 141-142           c
c  full description of algorithm in    cpc 21 (1981) 297-314           c
c  this version written up       in    cpc 27 (1982) 147-166           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       
      implicit real*8 (a-h,o-z)                                         
      dimension    fc(1),gc(1),fcp(1),gcp(1)                            
      logical      etane0,xlturn                                        
      common       /stee / paccq,nfp,npq,iexp,m1                        
      data zero,one,two,ten2,abort /0.0d0, 1.0d0, 2.0d0, 1.0d2, 2.0d5/  
      data half,tm30 / 0.5d0, 1.0d-30 /                                 
      data rt2dpi /0.79788 45608 02865 35587 98921 19868 76373 d0/      
c                                                                       
                        accur = 1.0d-16                                 
      mode  = 1                                                         
      if(mode1 .eq. 2 .or. mode1 .eq. 3 ) mode = mode1                  
      ifail = 0                                                         
      iexp  = 1                                                         
      npq   = 0                                                         
      eta   = eta1                                                      
      gjwkb = zero                                                      
      paccq = one                                                       
      if(kfn .ne. 0) eta = zero                                         
                 etane0  = eta .ne. zero                                
      acc   = accur * 10d0                                              
      acc4  = acc*ten2*ten2                                             
      acch  = dsqrt(acc)                                                
c ***    test range of xx, exit if.le.dsqrt(accur) or if negative       
c                                                                       
      if(xx .le. acch)                          go to 100               
      x     = xx                                                        
      xlm   = xlmin                                                     
      if(kfn .eq. 2)  xlm = xlm - half                                  
      if(xlm .le. -one .or. xlmax .lt. xlmin)   go to 105               
      e2mm1 = eta*eta + xlm*xlm + xlm                                   
      xlturn= x*(x - two*eta) .lt. xlm*xlm + xlm                        
      dell  = xlmax - xlmin + acc                                       
      if(dabs(dmod(dell,one)) .gt. acc) write(6,2040)xlmax,xlmin,dell   
      lxtra = idint(dell)                                               
      xll   = xlm + dfloat(lxtra)                                       
c ***       lxtra is number of additional lambda values to be computed  
c ***       xll  is max lambda value, or 0.5 smaller for j,y bessels    
c ***         determine starting array element (m1) from xlmin          
      m1  = max0(idint(xlmin + acc),0) + 1                              
      l1  = m1 + lxtra                                                  
c                                                                       
c ***    evaluate cf1  =  f   =  fprime(xl,eta,x)/f(xl,eta,x)           
c                                                                       
      xi  = one/x                                                       
      fcl = one                                                         
      pk  = xll + one                                                   
      px  = pk  + abort                                                 
      f   =  eta/pk + pk*xi                                             
         if(dabs(f).lt.tm30) f = tm30                                   
         d = zero                                                       
         c = f                                                          
c                                                                       
c ***   begin cf1 loop on pk = k = lambda + 1                           
c                                                                       
    4 pk1   = pk + one                                                  
        ek  = eta / pk                                                  
        rk2 = one + ek*ek                                               
        tk  = (pk + pk1)*(xi + ek/pk1)                                  
        d   =  tk - rk2 * d                                             
        c   =  tk - rk2 / c                                             
         if(dabs(c).lt.tm30) c = tm30                                   
         if(dabs(d).lt.tm30) d = tm30                                   
         d = one/d                                                      
         df = d * c                                                     
         f  = f * df                                                    
            if(d .lt. zero) fcl = - fcl                                 
         pk = pk1                                                       
                          if( pk .gt. px ) go to 110                    
      if(dabs(df-one) .ge. acc)             go to 4                     
                  nfp = pk - xll - 1                                    
      if(lxtra .eq. 0)                          go to 7                 
c                                                                       
c *** downward recurrence to lambda = xlm. array gc,if present,stores rl
c                                                                       
      fcl = fcl*tm30                                                    
      fpl = fcl*f                                                       
      if(mode .eq. 1) fcp(l1) = fpl                                     
                      fc (l1) = fcl                                     
      xl  = xll                                                         
      rl  = one                                                         
      el  = zero                                                        
      do 6  lp = 1,lxtra                                                
         if(etane0) el = eta/xl                                         
         if(etane0) rl = dsqrt(one + el*el)                             
         sl    =  el  + xl*xi                                           
         l     =  l1  - lp                                              
         fcl1  = (fcl *sl + fpl)/rl                                     
         fpl   =  fcl1*sl - fcl *rl                                     
         fcl   =  fcl1                                                  
         fc(l) =  fcl                                                   
         if(mode .eq. 1) fcp(l)  = fpl                                  
         if(mode .ne. 3 .and. etane0) gc(l+1) = rl                      
    6 xl = xl - one                                                     
      if(fcl .eq. zero) fcl = acc                                       
      f  = fpl/fcl                                                      
c ***    now we have reached lambda = xlmin = xlm                       
c ***    evaluate cf2 = p + i.q  again using steed's algorithm          
c ***    see text for compact complex code for sp cdc or non-ansi ibm   
c                                                                       
    7 if( xlturn ) call jwkb(x,eta,dmax1(xlm,zero),fjwkb,gjwkb,iexp)    
      if( iexp .gt. 1 .or. gjwkb .gt. one/(acch*ten2))  go to 9         
          xlturn = .false.                                              
      ta =  two*abort                                                   
      pk =  zero                                                        
      wi =  eta + eta                                                   
      p  =  zero                                                        
      q  =  one - eta*xi                                                
      ar = -e2mm1                                                       
      ai =  eta                                                         
      br =  two*(x - eta)                                               
      bi =  two                                                         
      dr =  br/(br*br + bi*bi)                                          
      di = -bi/(br*br + bi*bi)                                          
      dp = -xi*(ar*di + ai*dr)                                          
      dq =  xi*(ar*dr - ai*di)                                          
    8 p     = p  + dp                                                   
         q  = q  + dq                                                   
         pk = pk + two                                                  
         ar = ar + pk                                                   
         ai = ai + wi                                                   
         bi = bi + two                                                  
         d  = ar*dr - ai*di + br                                        
         di = ai*dr + ar*di + bi                                        
         c  = one/(d*d + di*di)                                         
         dr =  c*d                                                      
         di = -c*di                                                     
         a  = br*dr - bi*di - one                                       
         b  = bi*dr + br*di                                             
         c  = dp*a  - dq*b                                              
         dq = dp*b  + dq*a                                              
         dp = c                                                         
         if(pk .gt. ta)                         go to 120               
      if(dabs(dp)+dabs(dq).ge.(dabs(p)+dabs(q))*acc)   go to 8          
                      npq   = pk/two                                    
                      paccq = half*acc/dmin1(dabs(q),one)               
                      if(dabs(p) .gt. dabs(q)) paccq = paccq*dabs(p)    
c                                                                       
c *** solve for fcm = f at lambda = xlm,then find norm factor w=w/fcm   
c                                                                       
      gam = (f - p)/q                                                   
            if(q .le. acc4*dabs(p))             go to 130               
      w   = one/dsqrt((f - p)*gam + q)                                  
            go to 10                                                    
c *** arrive here if g(xlm) .gt. 10**6 or iexp .gt. 70 & xlturn = .true.
    9 w   = fjwkb                                                       
      gam = gjwkb*w                                                     
      p   = f                                                           
      q   = one                                                         
c                                                                       
c *** normalise for spherical or cylindrical bessel functions           
c                                                                       
   10                     alpha = zero                                  
          if(kfn  .eq. 1) alpha = xi                                    
          if(kfn  .eq. 2) alpha = xi*half                               
                          beta  = one                                   
          if(kfn  .eq. 1) beta  = xi                                    
          if(kfn  .eq. 2) beta  = dsqrt(xi)*rt2dpi                      
      fcm  = dsign(w,fcl)*beta                                          
           fc(m1)  = fcm                                                
                      if(mode .eq. 3)           go to 11                
           if(.not. xlturn)   gcl =  fcm*gam                            
           if(      xlturn)   gcl =  gjwkb*beta                         
           if( kfn .ne. 0 )   gcl = -gcl                                
           gc(m1)  = gcl                                                
           gpl =  gcl*(p - q/gam) - alpha*gcl                           
                      if(mode .eq. 2)           go to 11                
           gcp(m1) = gpl                                                
           fcp(m1) = fcm*(f - alpha)                                    
   11 if(lxtra .eq. 0 ) return                                          
c *** upward recurrence from gc(m1),gcp(m1)  stored value is rl         
c *** renormalise fc,fcp at each lambda and correct regular derivative  
c ***    xl   = xlm here  and rl = one , el = zero for bessels          
         w    = beta*w/dabs(fcl)                                        
         maxl = l1 - 1                                                  
      do 12 l = m1,maxl                                                 
                      if(mode .eq. 3)           go to 12                
                      xl = xl + one                                     
         if(etane0)   el = eta/xl                                       
         if(etane0)   rl = gc(l+1)                                      
                      sl = el + xl*xi                                   
         gcl1     = ((sl - alpha)*gcl - gpl)/rl                         
         gpl      =   rl*gcl -  (sl + alpha)*gcl1                       
         gcl      = gcl1                                                
         gc(l+1)  = gcl1                                                
                      if(mode .eq. 2)           go to 12                
         gcp(l+1) = gpl                                                 
         fcp(l+1) = w*(fcp(l+1) - alpha*fc(l+1))                        
   12 fc(l+1)     = w* fc(l+1)                                          
      return                                                            
c                                                                       
c ***    error messages                                                 
c                                                                       
  100 ifail = -1                                                        
      write(6,2000) xx,acch                                             
 2000 format(' for xx = ',1p,d12.3,' try small-x  solutions',           
     *' or x negative'/ ,' square root accuracy parameter =  ',d12.3/)  
      return                                                            
  105 ifail = -2                                                        
      write (6,2005) xlmax,xlmin,xlm                                    
 2005 format(/' problem with input order values:xlmax,xlmin,xlm = ',    
     *1p,3d15.6/)                                                       
      return                                                            
  110 ifail =  1                                                        
      write (6,2010) abort,f ,df,pk,px,acc                              
 2010 format(' cf1 has failed to converge after ',f10.0,' iterations',/ 
     *' f,df,pk,px,accur =  ',1p,5d12.3//)                              
      return                                                            
  120 ifail =  2                                                        
      write (6,2020) abort,p,q,dp,dq,acc                                
 2020 format(' cf2 has failed to converge after ',f7.0,' iterations',/  
     *' p,q,dp,dq,accur =  ',1p,4d17.7,d12.3//)                         
      return                                                            
  130 ifail =  3                                                        
      write (6,2030) p,q,acc,dell,lxtra,m1                              
 2030 format(' final q.le.dabs(p)*acc*10**4 , p,q,acc = ',1p,3d12.3,4x, 
     *' dell,lxtra,m1 = ',d12.3,2i5 /)                                  
      return                                                            
 2040 format(' xlmax - xlmin = dell not an integer ',1p,3d20.10/)       
      end                                                               
c                                                                       
      subroutine jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)                      
      real*8          xx,eta1,xl,fjwkb,gjwkb,dzero                      
c *** computes jwkb approximations to coulomb functions    for xl.ge. 0 
c *** as modified by biedenharn et al. phys rev 97 (1955) 542-554       
c *** calls dmax1,sqrt,alog,exp,atan2,float,int        barnett feb 1981 
      data   zero,half,one,six,ten/ 0.0e0, 0.5e0, 1.0e0, 6.0e0, 10.0e0 /
      data  dzero, rl35, aloge  /0.0d0, 35.0e0, 0.43429 45 e0 /         
      x     = xx                                                        
      eta   = eta1                                                      
      gh2   = x*(eta + eta - x)                                         
      xll1  = dmax1(xl*xl + xl,dzero)                                   
      if(gh2 + xll1 .le. zero) return                                   
       hll  = xll1 + six/rl35                                           
       hl   = sqrt(hll)                                                 
       sl   = eta/hl + hl/x                                             
       rl2  = one + eta*eta/hll                                         
       gh   = sqrt(gh2 + hll)/x                                         
       phi  = x*gh - half*( hl*alog((gh + sl)**2/rl2) - alog(gh) )      
          if(eta .ne. zero) phi = phi - eta*atan2(x*gh,x - eta)         
      phi10 = -phi*aloge                                                
      iexp  =  int(phi10)                                               
      if(iexp .gt. 70) gjwkb = ten**(phi10 - float(iexp))               
      if(iexp .le. 70) gjwkb = exp(-phi)                                
      if(iexp .le. 70) iexp  = 0                                        
      fjwkb = half/(gh*gjwkb)                                           
      return                                                            
      end                                                               
