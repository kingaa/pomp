      subroutine driverSSA(fprob,nvar,nevent,npar,nreps,ntimes,kflag,
     &     xstart,times,params,xout,e,v,d,nzero,izero,istate,ipar,
     &     ncovar,icovar,lcov,mcov,tcov,cov)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      external fprob
      ntreeh=1+int(log(nevent*1.0)/log(2.0)+1)
      do irep=1,nreps
         call SSA(fprob,irep,nvar,nevent,ntreeh,npar,nreps,
     &        ntimes,kflag,xstart,times,params,xout,e,v,d,
     &        nzero,izero,istate,ipar,ncovar,icovar,lcov,mcov,
     &        tcov,cov)
      enddo
      return
      end

      subroutine SSA(fprob,irep,nvar,nevent,ntreeh,npar,nreps,
     &     ntimes,kflag,xstart,times,params,xout,e,v,d,nzero,
     &     izero,istate,ipar,ncovar,icovar,lcov,mcov,tcov,cov)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension xstart(nvar,nreps), times(ntimes)
      dimension params(npar,nreps), par(npar)
      dimension xout(nvar,nreps,ntimes)
      dimension y(nvar),f(nevent,ntreeh),k(nevent)
      dimension e(nvar),v(nvar,nevent),d(nvar,nevent)
      dimension izero(nzero)
      dimension istate(nvar), ipar(npar), icovar(ncovar)
      dimension tcov(lcov), cov(lcov,mcov)
      dimension covars(mcov)
      external gillespie,kleap,fprob
      n=nvar
      m=nevent
c================
c Initialisation
c================
      t=times(1)
      tmax=times(ntimes)
      iflag=0
      icount=2
      do i=1,npar
         par(i)=params(i,irep)
      enddo
      do i=1,n
         y(i)=xstart(i,irep)
      enddo
c Set appropriate states to zero
      do i=1,nzero
         y(izero(i)+1)=0.0d0
      enddo
c Copy initial states into xout
      do i=1,n
         xout(i,irep,1)=y(i)
      enddo
c Initialize the covariate vector
      if(mcov.gt.0)then
         call tlook(lcov,mcov,tcov,cov,t,covars)
      endif
c========================================
c Initialise propensity functions & tree
c========================================
      do j=1,m
         f(j,1)=fprob(j,t,y,par,istate,ipar,icovar,mcov,covars)
      enddo
      jdum=m
      do itree=2,ntreeh
         jend=int((jdum+1)/2)
         do j=1,jend
            if(2*j.le.jdum)then
               f(j,itree)=f(2*j-1,itree-1)+f(2*j,itree-1)
            else
               f(j,itree)=f(2*j-1,itree-1)         
            endif
         enddo
         jdum=jend
      enddo
c=====
      do while(icount.le.ntimes)
         call rchkusr
         if(kflag.eq.0)then
            call gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,
     &           jevent,iflag,istate,ipar,ncovar,icovar,
     &           mcov,covars)
            if(iflag.eq.1)goto 100
         else                   !if(kflag.eq.1)then
c=================
c Determine kappa (most accurate but slowest method)
c=================
            dum=10d8
            do i=1,n
               dum=min(e(i)*y(i),dum)
               if(dum.le.1.0)goto 50
            enddo
 50         kappa=max(dum,1.0)
            if(kappa.eq.1)then
               call gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,
     &              jevent,iflag,istate,ipar,ncovar,icovar,
     &              mcov,covars)
               if(iflag.eq.1)goto 100
            else
               call kleap(fprob,kappa,t,f,y,v,d,par,n,m,ntreeh,npar,
     &              k,iflag,istate,ipar,ncovar,icovar,
     &              mcov,covars)
               if(iflag.eq.1)goto 100
            endif
c         else
c===============
c Determine tau (need to add code to avoid negative #s & determine tau)
c===============
c            tau=e(1)
c            call tauleap(fprob,tau,t,f,y,v,d,par,n,m,ntreeh,npar,
c     &           k,iflag)
c            if(iflag.eq.1)goto 100
         endif
c     
c Recording output at required time points
c     
         do while((icount.le.ntimes).and.(t.ge.times(icount)))
            do i=1,n
               xout(i,irep,icount)=y(i)
            enddo
c===============================
c     Set appropriate states to zero
c===============================
            do i=1,nzero
               y(izero(i)+1)=0.0d0
            enddo
            icount=icount+1
         enddo

         if((mcov.gt.0).and.(t.le.tmax))then
            call tlook(lcov,mcov,tcov,cov,t,covars)
         endif
c     
      enddo
c=====
 100  return
      end
      
      subroutine gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,jevent,
     &     iflag,istate,ipar,ncovar,icovar,mcov,cov)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension y(n),f(m,ntreeh),v(n,m),d(n,m),par(npar),ichangey(n)
      dimension istate(n),ipar(npar),icovar(ncovar),cov(mcov)
      external unifrnd,fprob
c=================================
c Generate uniform random numbers
c=================================
      p1=unifrnd()
      p2=unifrnd()
c=========================================
c Determine time interval and update time
c=========================================
      fsum=f(1,ntreeh)
      if(fsum.gt.0.0)then
         tstep=-log(p1)/fsum
         t=t+tstep
      else
         iflag=1
         goto 500
      endif
c=========================================
c Determine event, update pops & events
c=========================================
      jtree=1
      temp=p2*fsum
      do itree=ntreeh-1,1,-1
         if(itree.eq.1)then
            if(temp.lt.f(jtree,itree))then
               jevent=jtree
            else 
               jevent=jtree+1
            endif
         else
            if(temp.lt.f(jtree,itree))then
               jtree=2*jtree-1
            else
               temp=temp-f(jtree,itree)
               jtree=2*jtree+1
            endif
         endif
      enddo 
      do i=1,n
         ichangey(i)=0
      enddo
      do i = 1,n
         if(v(i,jevent).ne.0)then
            y(i) = y(i) + v(i,jevent)
            ichangey(i)=1
         endif
      enddo
c
c only updating events & tree entries that have changed
c
      do j=1,m
         do i=1,n
            if(ichangey(i).ne.0.and.d(i,j).ne.0)then
               fold=f(j,1)
               f(j,1)=fprob(j,t,y,par,istate,ipar,icovar,mcov,cov)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400     continue
      enddo
 500  return
      end

      subroutine kleap(fprob,kappa,t,f,y,v,d,par,n,m,ntreeh,npar,k,
     &     iflag,istate,ipar,ncovar,icovar,mcov,cov)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension y(n),f(m,ntreeh),p(m),v(n,m),d(n,m),par(npar)
      dimension k(m),ichangey(n)
      dimension istate(n),ipar(npar),icovar(ncovar),cov(mcov)
      external gammarnd,fprob
c=========================================
c Determine time interval and update time
c=========================================
      fsum=f(1,ntreeh)
      if(fsum.gt.0.0d0)then
         tstep=gammarnd(kappa,fsum)
         t=t+tstep
      else
         iflag=1
         goto 500
      endif
c=====================================================
c Determine frequency of events, update pops & events
c=====================================================
      do j=1,m
         p(j)=f(j,1)/fsum
      enddo
      call multinomrnd(kappa,p,m,k)
c
c some matrix-vector multiplication but only where necessary
c
      do i=1,n
         ichangey(i)=0
      enddo
      do j=1,m
         if(k(j).ne.0)then
            temp=k(j)
            do i = 1,n
               if(v(i,j).ne.0)then
                  y(i) = y(i) + temp*v(i,j)
                  ichangey(i)=1
               endif
            enddo
         endif
      enddo
c
c only updating events & tree entries that have changed
c
      do j=1,m
         do i=1,n
            if(ichangey(i).ne.0.and.d(i,j).ne.0)then
               fold=f(j,1)
               f(j,1)=fprob(j,t,y,par,istate,ipar,icovar,mcov,cov)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400     continue
      enddo
 500  return
      end
	  
      subroutine tauleap(fprob,tau,t,f,y,v,d,par,n,m,ntreeh,npar,k,
     &     iflag,istate,ipar,ncovar,icovar,mcov,cov)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension y(n),f(m,ntreeh),k(m),v(n,m),d(n,m),par(npar)
      dimension ichangey(n)
      dimension istate(n),ipar(npar),icovar(ncovar),cov(mcov)
      external poisrnd,fprob
c======================================================================
c Generate Poisson random variables (number of event firings in t+tau)
c======================================================================
      fsum=f(1,ntreeh)
      if(fsum.gt.0.0d0)then
         do j=1,m
            k(j)=poisrnd(f(j,1)*tau)
         enddo
      else
         iflag=1
         goto 500
      endif
c==========
c Update t
c==========
      t=t+tau
c================================================
c Compute changes in population numbers & events
c================================================
c
c some matrix-vector multiplication but only where necessary
c
      do i=1,n
         ichangey(i)=0
      enddo
      do j=1,m
         if(k(j).ne.0)then
            temp=k(j)
            do i = 1,n
               if(v(i,j).ne.0)then
                  y(i) = y(i) + temp*v(i,j)
                  ichangey(i)=1
               endif
            enddo
         endif
      enddo
c
c only updating events & tree entries that have changed
c
      do j=1,m
         do i=1,n
            if(ichangey(i).ne.0.and.d(i,j).ne.0)then
               fold=f(j,1)
               f(j,1)=fprob(j,t,y,par,istate,ipar,icovar,mcov,cov)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400     continue
      enddo
 500  return
      end
