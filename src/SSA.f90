subroutine driverSSA(fprob,nvar,nevent,npar,nreps,ntimes,kflag,&
     xstart,times,params,xout,e,v,d,ndeps,ideps,&
     nzero,izero,istate,ipar,ncovar,&
     icovar,lcov,mcov,tcov,cov,iflag)
  implicit integer (i-n)
  implicit double precision (a-h,o-z)
  dimension xstart(nvar,nreps), times(ntimes)
  dimension params(npar,nreps)
  dimension xout(nvar,nreps,ntimes)
  dimension e(nvar),v(nvar,nevent),d(nvar,nevent),ideps(ndeps)
  dimension izero(nzero)
  dimension istate(nvar), ipar(npar), icovar(ncovar)
  dimension tcov(lcov), cov(lcov,mcov)
  external fprob
  ntreeh=1+int(log(nevent*1.0)/log(2.0)+1)
  if(iflag.ne.0)return
  do irep=1,nreps
     call SSA(fprob,irep,nvar,nevent,ntreeh,npar,nreps,&
          ntimes,kflag,xstart,times,params,xout,e,v,d,ndeps,ideps,&
          nzero,izero,istate,ipar,ncovar,icovar,lcov,mcov,&
          tcov,cov,iflag)
     if(iflag.eq.2)return
  enddo
  return
end subroutine driverSSA

subroutine SSA(fprob,irep,nvar,nevent,ntreeh,npar,nreps,ntimes,&
     kflag,xstart,times,params,xout,e,v,d,ndeps,ideps,nzero,&
     izero,istate,ipar,ncovar,icovar,lcov,mcov,tcov,cov,iflag)
  implicit integer (i-n)
  implicit double precision (a-h,o-z)
  dimension xstart(nvar,nreps), times(ntimes)
  dimension params(npar,nreps), par(npar)
  dimension xout(nvar,nreps,ntimes)
  dimension y(nvar),f(nevent,ntreeh),k(nevent)
  dimension e(nvar),v(nvar,nevent),d(nvar,nevent),ideps(ndeps)
  dimension izero(nzero)
  dimension istate(nvar), ipar(npar), icovar(ncovar)
  dimension tcov(lcov), cov(lcov,mcov)
  dimension covars(mcov)
  external gillespie,kleap,fprob
  n=nvar
  m=nevent
  !================
  ! Initialisation
  !================
  t=times(1)
  tmax=times(ntimes)
  icount=2
  do i=1,npar
     par(i)=params(i,irep)
  enddo
  do i=1,n
     y(i)=xstart(i,irep)
  enddo
  ! Set appropriate states to zero
  do i=1,nzero
     y(izero(i)+1)=0.0d0
  enddo
  ! Copy initial states into xout
  do i=1,n
     xout(i,irep,1)=y(i)
  enddo
  ! Initialize the covariate vector
  if(mcov.gt.0)then
     call tlook(lcov,mcov,tcov,cov,t,covars)
  endif
  !========================================
  ! Initialise propensity functions & tree
  !========================================
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
  !=====
  do while(icount.le.ntimes)
     call rchkusr
     if(kflag.eq.0)then
        call gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,&
             jevent,iflag,istate,ipar,ncovar,icovar,&
             mcov,covars)
     else
        !=================
        ! Determine kappa (most accurate but slowest method)
        !=================
        dum=10d8
        do i=1,ndeps
           dum=min(e(ideps(i))*y(ideps(i)),dum)
           if(dum.le.1.0)goto 50
        enddo
50      kappa=int(max(dum,1.0d0))
        if(kappa.eq.1)then
           call gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,&
                jevent,iflag,istate,ipar,ncovar,icovar,&
                mcov,covars)
        else
           call kleap(fprob,kappa,t,f,y,v,d,par,n,m,ntreeh,npar,&
                k,iflag,istate,ipar,ncovar,icovar,&
                mcov,covars)
        endif
     endif
     if(iflag.eq.2)goto 100
     !
     ! Recording output at required time points
     !
     do while ((icount.le.ntimes).and.((iflag.eq.1).or.(t.ge.times(icount))))
        do i=1,n
           xout(i,irep,icount)=y(i)
        enddo
        !===============================
        !     Set appropriate states to zero
        !===============================
        do i=1,nzero
           y(izero(i)+1)=0.0d0
        enddo
        if (iflag.eq.1)then
           t = times(icount)
        endif
        icount=icount+1
        if (icount.gt.ntimes) exit
     enddo

     if((mcov.gt.0).and.(t.le.tmax))then
        call tlook(lcov,mcov,tcov,cov,t,covars)
     endif
     !
  enddo
  !=====
100 return
end subroutine SSA

subroutine gillespie(fprob,t,f,y,v,d,par,n,m,ntreeh,npar,jevent,&
     iflag,istate,ipar,ncovar,icovar,mcov,cov)
  implicit integer (i-n)
  implicit double precision (a-h,o-z)
  dimension y(n),f(m,ntreeh),v(n,m),d(n,m),par(npar),ichangey(n)
  dimension istate(n),ipar(npar),icovar(ncovar),cov(mcov)
  external unifrnd,fprob
  !=================================
  ! Generate uniform random numbers
  !=================================
  p1=unifrnd()
  p2=unifrnd()
  !=========================================
  ! Determine time interval and update time
  !=========================================
  fsum=f(1,ntreeh)
  if(fsum>0.0d0)then
     iflag=0
     tstep=-log(p1)/fsum
     t=t+tstep
  elseif(fsum<0.0d0)then
     iflag=2
     goto 500
  else
     iflag=1
     goto 500
  endif
  !=========================================
  ! Determine event, update pops & events
  !=========================================
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
  !
  ! only updating events & tree entries that have changed
  !
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
400  continue
  enddo
500 return
end subroutine gillespie

subroutine kleap(fprob,kappa,t,f,y,v,d,par,n,m,ntreeh,npar,k,&
     iflag,istate,ipar,ncovar,icovar,mcov,cov)
  implicit integer (i-n)
  implicit double precision (a-h,o-z)
  dimension y(n),f(m,ntreeh),p(m),v(n,m),d(n,m),par(npar)
  dimension k(m),ichangey(n)
  dimension istate(n),ipar(npar),icovar(ncovar),cov(mcov)
  external gammarnd,multinomrnd,fprob
  !=========================================
  ! Determine time interval and update time
  !=========================================
  fsum=f(1,ntreeh)
  if(fsum>0.0d0)then
     iflag=0
     tstep=gammarnd(kappa,fsum)
     t=t+tstep
  elseif(fsum<0.0d0)then
     iflag=2
     goto 500
  else
     iflag=1
     goto 500
  endif
  !=====================================================
  ! Determine frequency of events, update pops & events
  !=====================================================
  do j=1,m
     p(j)=f(j,1)/fsum
  enddo
  call multinomrnd(kappa,p,m,k)
  !
  ! some matrix-vector multiplication but only where necessary
  !
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
  !
  ! only updating events & tree entries that have changed
  !
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
400  continue
  enddo
500 return
end subroutine kleap
