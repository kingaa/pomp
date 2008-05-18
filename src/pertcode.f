      subroutine pertcode(xstart,times,params,ngams,xout,nvar, 
     &     npar,nreps,ntimes)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension xstart(nvar,nreps), times(ntimes), params(npar,nreps)
      dimension xout(nvar+2,nreps,ntimes)
      dimension ngams(6)
      do irep=1,nreps
         npE1=ngams(1)
         npI1=ngams(2)
         npR=ngams(3)
         npE2=ngams(4)
         npI2=ngams(5)
         npV=ngams(6)
         nsumgam=npE1+npI1+npR+npE2+npI2+npV
c         if(nvar.ne.nsumgam+2)then
c            print *,'error: mismatch in state variable dimension'
c            stop
c         endif
         nstore=nvar+3
         nevent=2*nsumgam+6
         ntreeh=1+int(log(nevent*1.0)/log(2.0)+1)
         call pertussis(irep,npE1,npI1,npR,npE2,npI2,npV,nvar,nstore,
     &        nevent,ntreeh,npar,nreps,ntimes,xstart,times,params,xout)
      enddo
      return
      end

      subroutine pertussis(irep,npE1,npI1,npR,npE2,npI2,npV,nvar,nstore,
     &        nevent,ntreeh,npar,nreps,ntimes,xstart,times,params,xout)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension xstart(nvar,nreps), times(ntimes), params(npar,nreps)
      dimension xout(nvar+2,nreps,ntimes)
      dimension y(nstore),yobs(2),e(nstore),f(nevent,ntreeh),k(nevent)
      dimension v(nstore,nevent),d(nstore,nevent)
      common /aliases/iS1,iS2,iI1,iI2,iN
      common /seasonality/beta11,b1
      common /transmission/beta12,beta21,beta22
      common /infection/sigma1,sigma2,gamma1,gamma2
      common /distn_params1/nE1,nI1,nR,nE2,nI2,nV
      common /distn_params2/nsumnI1,nsumnR,nsumnE2,nsumnI2,nsumnV
      common /immigration/rho1,rho2
      common /immunity/alphan,alphav
      common /vitals/rnu,rmu,vp
      external gillespie,kleap
      n=nstore
      m=nevent
      nE1=npE1
      nI1=npI1
      nR=npR
      nE2=npE2
      nI2=npI2
      nV=npV
c==================
c Variable aliases
c==================
      iS1=1
      iS2=iS1+nE1+nI1+nR+1
      iI1=n-2
      iI2=n-1
      iN=n
      nsumnI1=nE1+nI1
      nsumnR=nsumnI1+nR
      nsumnE2=nsumnR+nE2
      nsumnI2=nsumnE2+nI2
      nsumnV=nsumnI2+nV
c=================================
c Read and store parameter values
c=================================
      sigma1 = params(1,irep)
      sigma2 = params(2,irep)
c recovery rate
      gamma1 = params(3,irep)
      gamma2 = params(4,irep)
c vital rate 
      rnu = params(5,irep)
      rmu = params(6,irep)
c primary transmission
      beta11 = params(7,irep)
      b1 = params(8,irep)
c contact structure
      beta12 = params(9,irep)
c secondary transmission
      beta21 = params(10,irep)
      beta22 = params(11,irep)
c immunity
      alphan = params(12,irep)
      alphav = params(13,irep)
c vaccination
      vp = params(14,irep)
c background influx of infections per year
      rho1 = params(15,irep)
      rho2 = params(16,irep)
c initial population size
      rinitN = params(17,irep)
c simulation params
      tmax = times(ntimes)
      epsilon=params(18,irep)
      t0=params(19,irep)
c
      do i=1,n
	e(i)=epsilon	
      enddo	
      e(iS1)=epsilon/2
      e(iS2)=epsilon/2
      e(iI1)=epsilon/2
      e(iI2)=epsilon/2
c
      do i=1,n
         do j=1,m
            v(i,j)=0
            if(j.eq.3.and.b1.gt.0.0)then
               d(i,j)=1 ! always update 1st transmission event if seasonal
            else
               d(i,j)=0
            endif
         enddo
      enddo
c
      v(iS1,1)=1
      v(iS1,(/2,3/))=-1
      d(iS1,(/2,3/))=1
      v(iN,2)=-1
c
      do i=2,nE1+1
         j=i+2
         if(i.eq.2)then
            v(i,j-1)=1
         else
            v(i,j+nE1-1)=1
         endif
         v(i,(/j,j+nE1/))=-1
         d(i,(/j,j+nE1/))=1
         v(iN,j)=-1
      enddo
      do i=nE1+2,nsumnI1+1
         j=i+2+nE1
         if(i.eq.nE1+2)then
            v(i,j-1)=1
         else
            v(i,j+nI1-1)=1
         endif
         v(i,(/j,j+nI1/))=-1
         d(i,(/j,j+nI1/))=1
         v(iI1,j)=-1
         v(iN,j)=-1
      enddo
      do i=nsumnI1+2,nsumnR+1
         j=i+2+nsumnI1
         if(i.eq.nsumnI1+2)then
            v(i,(/j-1,2*(nsumnI2+2)+1,m/))=1
            v(iI1,(/j+nR-2/))=-1
            v(iI2,(/2*(nsumnI2+2)+1/))=-1
         else
            v(i,j+nR-1)=1
         endif
         v(i,(/j,j+nR/))=-1
         d(i,(/j,j+nR/))=1
         v(iN,j)=-1
      enddo
      v(iS2,2*(nsumnR+1)+1)=1
      v(iS2,(/2*(nsumnR+2),2*(nsumnR+2)+1/))=-1
      d(iS2,(/2*(nsumnR+2),2*(nsumnR+2)+1/))=1
      v(iN,2*(nsumnR+2))=-1
      do i=iS2+1,iS2+nE2
         j=i+3+nsumnR
         if(i.eq.iS2+1)then
            v(i,j-1)=1
         else
            v(i,j+nE2-1)=1
         endif
         v(i,(/j,j+nE2/))=-1
         d(i,(/j,j+nE2/))=1
         v(iN,j)=-1
      enddo
      do i=iS2+nE2+1,iS2+nE2+nI2
         j=i+3+nsumnE2
         if(i.eq.iS2+nE2+1)then
            v(i,j-1)=1
         else
            v(i,j+nI2-1)=1
         endif
         v(i,(/j,j+nI2/))=-1
         d(i,(/j,j+nI2/))=1
         v(iI2,j)=-1
         v(iN,j)=-1
      enddo
      do i=iS2+nE2+nI2+1,iS2+nE2+nI2+nV
         j=i+4+nsumnI2
         if(i.eq.iS2+nE2+nI2+1)then
            v(i,j-1)=1
         else
            v(i,j+nV-1)=1
         endif
         v(i,(/j,j+nV/))=-1
         d(i,(/j,j+nV/))=1
         v(iN,j)=-1
      enddo
      v(iI1,3+2*nE1)=1
      d(iI1,(/3,2*(nsumnR+2)+1/))=1
      v(iI2,3+2*(nsumnE2+1))=1
      d(iI2,(/3,2*(nsumnR+2)+1/))=1
      v(iN,(/1,2*(nsumnI2+3)/))=1
      d(iN,(/1,3,2*(nsumnR+2)+1,2*(nsumnI2+3)/))=1
c=====
c      iseed1=500 !time() ! seed for random number generator
c      iseed2=7395+iseed1
c      call setall(iseed1,iseed2)
c================
c Initialization
c================
      t=t0 ! allow for burn-in if t0<0
      iflag=0
      icount=1
      crpnew=0.0
      crsnew=0.0
      do i=1,nvar
         y(i)=xstart(i,irep)
      enddo
      ysumI1=0.0
      do i=iS1+nE1+1,iS1+nE1+nI1
         ysumI1=ysumI1+y(i)
      enddo
      ysumI2=0.0
      do i=iS2+nE2+1,iS2+nE2+nI2
         ysumI2=ysumI2+y(i)
      enddo
      y(iI1)=ysumI1
      y(iI2)=ysumI2
      y(iN)=rinitN
c========================================
c Initialise propensity functions & tree
c========================================
      do j=1,m
         f(j,1)=fprob(j,n,t,y)
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
      do while(t.lt.tmax)
c=================
c Determine kappa
c=================
         dum=e(1)*y(iS1)
         do i=2,n
            dum=min(e(i)*y(i),dum)
         enddo
         kappa=max(dum,1.0)
         if(kappa.eq.1)then
            call gillespie(t,f,y,v,d,n,m,ntreeh,jevent,iflag)
            if(iflag.eq.1)goto 100
            if(jevent.eq.3+2*nsumnI1)crpnew=crpnew+1
            if(jevent.eq.5+2*nsumnI2)crsnew=crsnew+1
         else
            call kleap(kappa,t,f,y,v,d,n,m,ntreeh,k,iflag)
            if(iflag.eq.1)goto 100
            crpnew=crpnew+k(3+2*nsumnI1)
            crsnew=crsnew+k(5+2*nsumnI2)
         endif
c
c Recording output at required time points
c
c Because t0<times(1), resetting case report accumulation in interval
c prior to first required output
c
         if(t.lt.times(1).and.t.ge.(2*times(1)-times(2)))then
            crpnew=0.0d0
            crsnew=0.0d0
         endif           
         if(t.ge.times(icount))then
            do i=ivar,nvar
               xout(ivar,irep,icount)=y(i)
            enddo
            xout(nvar+1,irep,icount)=crpnew
            xout(nvar+2,irep,icount)=crsnew
            icount=icount+1
            crpnew=0.0d0
            crsnew=0.0d0
         endif
c
      enddo
c=====
 100  return
      end
      
      function fprob(j,n,t,y)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension y(n)
      common /aliases/iS1,iS2,iI1,iI2,iN
      common /transmission/beta12,beta21,beta22
      common /infection/sigma1,sigma2,gamma1,gamma2
      common /distn_params1/nE1,nI1,nR,nE2,nI2,nV
      common /distn_params2/nsumnI1,nsumnR,nsumnE2,nsumnI2,nsumnV
      common /immigration/rho1,rho2
      common /immunity/alphan,alphav
      common /vitals/rnu,rmu,vp
c
      yFI1 = (beta(t)*y(iI1)+beta12*y(iI2)+rho1)/y(iN)
      yFI2 = (beta21*y(iI1)+beta22*y(iI2)+rho2)/y(iN)
c
      if(j.eq.1)fprob = (1.0d0-vp)*rnu*y(iN) ! Birth
      if(j.eq.2)fprob = rmu*y(iS1)         ! Susceptible death
      if(j.eq.3)fprob = yFI1*y(iS1)          ! Infection
c     
      if(j.gt.3.and.j.le.(3+nE1))fprob = rmu*y(iS1+j-3) ! Exposed death
      if(j.gt.(3+nE1).and.j.le.(3+2*nE1))
     &     fprob = sigma1*nE1*y(iS1+j-(3+nE1)) ! Movement into next class
c
      if(j.gt.(3+2*nE1).and.j.le.(3+nE1+nsumnI1))
     &     fprob = rmu*y(iS1+j-(3+nE1))      ! Infected death
      if(j.gt.(3+nE1+nsumnI1).and.j.le.(3+2*nsumnI1))
     &     fprob = gamma1*nI1*y(iS1+j-(3+nsumnI1))   ! Movement into next class
c
      if(j.gt.(3+2*nsumnI1).and.j.le.(3+nsumnI1+nsumnR))
     &     fprob = rmu*y(iS1+j-(3+nsumnI1))      ! Recovered death
      if(j.gt.(3+nsumnI1+nsumnR).and.j.le.(3+2*nsumnR))
     &     fprob =  alphan*nR*y(iS1+j-(3+nsumnR))   ! Movement into next class
c
      if(j.eq.4+2*nsumnR)fprob = rmu*y(iS2)        ! Susceptible death
      if(j.eq.5+2*nsumnR)fprob = yFI2*y(iS2)         ! Infection
c     
      if(j.gt.5+2*nsumnR.and.j.le.5+nsumnR+nsumnE2)
     &     fprob = rmu*y(iS2+j-(5+2*nsumnR))     ! Exposed death
      if(j.gt.5+nsumnR+nsumnE2.and.j.le.5+2*nsumnE2)
     &     fprob = sigma2*nE2*y(iS2+j-(5+nsumnR+nsumnE2))   ! Movement into next class
c     
      if(j.gt.5+2*nsumnE2.and.j.le.5+nsumnE2+nsumnI2)
     &     fprob = rmu*y(iS2+j-(5+nsumnR+nsumnE2))     ! Infected death
      if(j.gt.5+nsumnE2+nsumnI2.and.j.le.5+2*nsumnI2)
     &     fprob = gamma2*nI2*y(iS2+j-(5+nsumnR+nsumnI2))  ! Movement into next class
c     
      if(j.eq.6+2*nsumnI2)fprob = vp*rnu*y(iN)      ! Vaccinated Birth
      if(j.gt.6+2*nsumnI2.and.j.le.6+nsumnI2+nsumnV)
     &     fprob = rmu*y(iS2+j-(6+nsumnR+nsumnI2))      ! Vaccinated death
      if(j.gt.6+nsumnI2+nsumnV.and.j.le.6+2*nsumnV)
     &     fprob = alphav*nV*y(iS2+j-(6+nsumnR+nsumnV))   ! Movement into next class
c=====
      return
      end

      function beta(t)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      common /seasonality/beta11,b1
      nyear=int(t)
      tdays=(t-nyear)*365.0d0
      if((tdays.ge.7.and.tdays.lt.100).or.(tdays.ge.116.and.
     &     tdays.lt.200).or.(tdays.ge.252.and.tdays.lt.300).or.
     &     (tdays.ge.308.and.tdays.lt.356))then
         iterm=1
      else
         iterm=-1
      endif
      b0 = beta11/(1+b1*181.0d0/365.0d0)
      beta=b0*(1.0d0+iterm*b1)
      return 
      end

