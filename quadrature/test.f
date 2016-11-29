      program main
      double precision a,b,result,abserr,resabs,resasc
      double precision book1,book3,book11,book15,book16
      double precision book454,book455,book457,book458,book459
      double precision myfn1,myfn2,myfn3,myfn4
      double precision alpha,beta
      double precision alist(1000),blist(1000),rlist(1000)
      double precision elist(1000),pts(1000)
      double precision rslst(1000),erlst(1000)
      double precision points(4)
      double precision chebmo(1000,25)
      integer iord(1000)
      integer nnlog(1000)
      integer ndin(1000)
      integer level(1000)
      integer ierlst(1000)
      integer maxp1,momcom
      integer inf
      integer integr
      common /ALPHA/alpha
      external book1,book3,book11,book15,book16
      external book454,book455,book457,book458,book459
      external myfn1,myfn2,myfn3,myfn4
      call gsl_ieee_env_setup

      a = 0.0
      b = 1.0

c      alpha = 2.6
c      print *,'alpha = ',alpha
c      call dqk15(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk21(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk31(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk41(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk51(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk61(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc

c      alpha = -0.9
c      print *,'alpha = ',alpha
c      call dqk15(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk21(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk31(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk41(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk51(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk61(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc


c     alpha = 1.3
      a = 0.3
      b = 2.71
c     print *,'OSCILLATINg alpha = ',alpha
c     call dqk15(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk21(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk31(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk41(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk51(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk61(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc



c      ier = 0
c      neval = 0
c      alpha = 2.6
c      epsabs = 1.0e-1
c      epsrel = 0.0
c      a = 0.0
c      b = 1.0
c      print *,'QNG book1'
c      print *,'alpha = ',alpha
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      epsabs = 0.0
c      epsrel = 1e-9
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      epsabs = 0.0
c      epsrel = 1e-13
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      alpha = -0.9
c      epsabs = 0.0
c      epsrel = 1e-3
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier

c     alpha = 1.3
c      a = 0.3
c      b = 2.71
c      epsabs = 0.0
c      epsrel = 1e-12
c      call dqng(book3,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier

c      alpha = 2.6
c      a = 0.0
c      b = 1.0
c      epsabs = 0.0
c      epsrel = 1d-10
c      key = 1
c      limit = 1000
c      print *, 'DQAGE'
c      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do 10 i=1,10
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c 10   continue

c      alpha = 2.6
c      a = 0.0
c      b = 1.0
c      epsabs = 1d-14
c      epsrel = 0
c      key = 2
c      limit = 1000
c      print *, 'DQAGE'
c      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do 11 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c 11   continue
      
c     alpha = 1.3
c     a = 0.3
c     b = 2.71
c     epsabs = 1d-14
c     epsrel = 0
c     key = 3
c     limit = 1000
c     print *, 'DQAGE oscill'
c     call dqage(book3,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 12 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c12   continue

c     alpha = 2.0
c     a = -1.0
c     b =  1.0
c     epsabs = 1d-14
c     epsrel = 0
c     key = 5
c     limit = 1000
c     print *, 'DQAGE sing'
c     call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 13 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c13   continue


c     alpha = 1.0
c     a = -1.0
c     b =  1.0
c     epsabs = 1d-7
c     epsrel = 0
c     key = 6
c     limit = 3
c     print *, 'DQAGE sing'
c     call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 14 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c14   continue


c     alpha = 2.6
c     a = 0.0
c     b = 1.0
c     epsabs = 0.0
c     epsrel = 1d-10
c     limit = 1000
c     print *, 'DQAGSE'
c     call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 15 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c15   continue


c     alpha = 2.6
c     a = 0.0
c     b = 1.0
c     epsabs = 1d-14
c     epsrel = 0.0
c     limit = 1000
c     print *, 'DQAGSE abs'
c     call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 16 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c16   continue

c      alpha = 2.0
c      a = 1.0
c      b = 1000.0
c      epsabs = 1d-7
c      epsrel = 0.0
c      limit = 1000
c      print *, 'DQAGSE abs'
c      call dqagse(book11,a,b,epsabs,epsrel,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,10
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c      alpha = 2.0
c      a = 0.0
c      inf = 1
c      epsabs = 0.0
c      epsrel = 1.0d-3
c      limit = 1000
c      print *, 'DQAGI abs'
c      call dqagie(book455,a,inf,epsabs,epsrel,limit,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,20
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c      alpha = 5.0
c      a = 0.0
c      inf = 1
c      epsabs = 0.0
c      epsrel = 1.0d-7
c      limit = 1000
c      print *, 'DQAGI abs'
c      call dqagie(book15,a,inf,epsabs,epsrel,limit,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,20
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c      alpha = 1.0
c      a = 99.9
c      inf = 1
c      epsabs = 1.0d-7
c      epsrel = 0.0
c      limit = 1000
c      print *, 'DQAGI abs'
c      call dqagie(book16,a,inf,epsabs,epsrel,limit,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,20
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c      alpha = 1.0
c      a = 0.0
c      inf = 2
c      epsabs = 1.0d-7
c      epsrel = 0.0
c      limit = 1000
c      print *, 'DQAGI abs'
c      call dqagie(myfn1,a,inf,epsabs,epsrel,limit,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,20
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c      alpha = 1.0
c      b = 1.0
c      inf = -1
c      epsabs = 1.0d-7
c      epsrel = 0.0
c      limit = 1000
c      print *, 'DQAGI abs'
c      call dqagie(myfn2,b,inf,epsabs,epsrel,limit,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,20
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo
      
c     alpha = 1.0
c     a = 0.0
c     b = 3.0
c     npts2 = 4
c     points(1) = 1.0
c     points(2) = sqrt(2.0d0)
c     epsabs = 0.0
c     epsrel = 1.0d-3
c     limit = 1000
c     print *, 'DQAGP'
c     call dqagpe(book454,a,b,npts2,points,epsabs,epsrel,limit,
c    $     result, abserr, neval,ier,alist,blist,rlist,elist,
c    $     pts,iord,level,ndin,last)
c     write(6,3) result, abserr, neval, ier, last
c     do i=1,25
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c     enddo

c      alpha = 1.0
c      a = -1.0
c      b = 5.0
c      c = 0.0
c      epsabs = 0.0
c      epsrel = 1.0d-3
c      limit = 1000
c      print *, 'DQAGP'
c      call dqawce(book459,a,b,c,epsabs,epsrel,limit,
c     $     result, abserr, neval,ier,alist,blist,rlist,elist,
c     $     iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,25
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

c     alpha = 0.0
c     beta = 0.0
c     integr = 2
c     a = 0.0
c     b = 1.0
c     epsabs = 0.0
c     epsrel = 1.0d-7
c     limit = 1000
c     print *, 'DQAGP'
c     call dqawse(book458,a,b,alpha,beta,integr,epsabs,epsrel,limit,
c    $     result, abserr, neval,ier,alist,blist,rlist,elist,
c    $     iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do i=1,25
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c     enddo

c     alpha = -0.5
c     beta = -0.3
c     integr = 4
c     a = 0.0
c     b = 1.0
c     epsabs = 0.0
c     epsrel = 1.0d-7
c     limit = 1000
c     print *, 'DQAGP'
c     call dqawse(book458,a,b,alpha,beta,integr,epsabs,epsrel,limit,
c    $     result, abserr, neval,ier,alist,blist,rlist,elist,
c    $     iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c     enddo

c     alpha = -0.5
c     beta = -0.3
c     integr = 4
c     a = 0.0
c     b = 1.0
c     epsabs = 0.0
c     epsrel = 1.0d-7
c     limit = 1000
c     print *, 'DQAGP'
c     call dqawse(book458,a,b,alpha,beta,integr,epsabs,epsrel,limit,
c    $     result, abserr, neval,ier,alist,blist,rlist,elist,
c    $     iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c     enddo

c     a = 0.0
c     b = 1.0
c     omega = 10.0 * 3.14159265358979323846
c     epsabs = 0.0
c     epsrel = 1.0d-7
c     limit = 1000
c     integr = 2
c     icall = 1
c     maxp1 = 1000
c     momcom = 0
c     print *, 'DQAGP'
c     call dqawoe(myfn4,a,b,omega,integr,epsabs,epsrel,limit,
c    $     icall, maxp1,
c    $     result, abserr, neval,ier,last,alist,blist,rlist,elist,
c    $     iord, nnlog, momcom, chebmo)
c     write(6,3) result, abserr, neval, ier, last
c     do i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c     enddo

      a = 0.0
      b = 1.0
      omega = 3.14159265358979323846 / 2.0
      epsabs = 1.0d-7
      epsrel = 0.0
      limit = 1000
      limlst = 1000
      integr = 1
      icall = 1
      maxp1 = 1000
      momcom = 0
      print *, 'DQAGP'
      call dqawfe(book457,a,omega,integr,epsabs,limlst, limit,
     $     maxp1,
     $     result, abserr, neval,ier,rslst, erlst, ierlst, lst,
     $     alist,blist,rlist,elist,
     $     iord, nnlog, chebmo)
      write(6,3) result, abserr, neval, ier, last
      do i=1,20
         write(6,4) i,rslst(i),erlst(i)
      enddo


 1    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'double exp_resabs =',1pe25.18, ';',/
     $     'double exp_resasc =',1pe25.18, ';')
 2    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'double exp_neval  =',I8, ';',/
     $     'double exp_ier    =',I8, ';')
 3    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'int    exp_neval  =',I8, ';',/
     $     'int    exp_ier    =',I8, ';',/
     $     'int    last       =',I8, ';')
 4    format('i=',i4,' a=',1pe25.18,' b=',1pe25.18,
     $     ' r=',1pe25.18,' e=',1pe25.18, ' iord=',i4)
      end




      double precision function book1(x)
      double precision alpha,x
      common /ALPHA/alpha
      book1=x**alpha*log(1.0/x)
      end

      double precision function book3(x)
      double precision alpha,x
      common /ALPHA/alpha
      book3=cos((2**alpha)*sin(x))
      end

      double precision function book11(x)
      double precision alpha,x
      common /ALPHA/alpha
      book11=(log(1/x))**(alpha-1.0)
      end

      double precision function book15(x)
      double precision alpha,x
      common /ALPHA/alpha
      book15=(x**2.0)*exp(-x*(2.0**-alpha))
      end

      double precision function book16(x)
      double precision alpha,x
      common /ALPHA/alpha
      book16=(x**(alpha-1.0))/((1.0+10.0*x)**2.0)
      end

      double precision function book454(x)
      double precision alpha,x
      common /ALPHA/alpha
      book454=(x**3.0) * log(abs((x**2.0 - 1.0)*(x**2.0 - 2.0)))
      write(6,6661) x, book454
 6661 format("FF x = ", 1pe25.18, " book454 = ", 1pe25.18)
      end

      double precision function book455(x)
      double precision alpha,x
      common /ALPHA/alpha
      book455=log(x)/(1.0 + 100.0*x*x)
      write(6,6661) x, book455
 6661 format("FF x = ", 1pe25.18, " book455 = ", 1pe25.18)
      end

      double precision function book457(x)
      double precision alpha,x
      common /ALPHA/alpha
      book457 = 0
      if (x.gt.0) book457= 1.0/sqrt(x)
      write(6,6661) x, book457
 6661 format("FF x = ", 1pe25.18, " book457 = ", 1pe25.18)
      end


      double precision function book458(x)
      double precision alpha,x
      common /ALPHA/alpha
      book458=1.0/(1.0 + log(x)**2)**2
      write(6,6661) x, book458
 6661 format("FF x = ", 1pe25.18, " book458 = ", 1pe25.18)
      end


      double precision function book459(x)
      double precision alpha,x
      common /ALPHA/alpha
      book459=1.0/(5.0*x**3 + 6.0)
      write(6,6661) x, book459
 6661 format("FF x = ", 1pe25.18, " book459 = ", 1pe25.18)
      end


      double precision function myfn1(x)
      double precision alpha,x
      common /ALPHA/alpha
      myfn1=exp(-x - x**2.0)
      end

      double precision function myfn2(x)
      double precision alpha,x
      common /ALPHA/alpha
      myfn2=exp(alpha*x)
      end

      double precision function myfn3(x)
      double precision alpha,x
      common /ALPHA/alpha
      myfn3=exp(-x)
      write(6,6661) x, myfn3
 6661 format("FF x = ", 1pe25.18, " myfn3 = ", 1pe25.18)
      end

      double precision function myfn4(x)
      double precision alpha,x
      common /ALPHA/alpha
      myfn4=0
      if (x.gt.0.0) myfn4=log(x)
      write(6,6661) x, myfn4
 6661 format("FF x = ", 1pe25.18, " myfn4 = ", 1pe25.18)
      end





