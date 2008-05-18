      subroutine gillespie(t,f,y,v,d,n,m,ntreeh,jevent,iflag)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      real ranf
      dimension y(n),f(m,ntreeh),v(n,m),d(n,m),ichangey(n)
      external ranf,fprob
c=================================
c Generate uniform random numbers
c=================================
      p1=ranf()
      p2=ranf()
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
               f(j,1)=fprob(j,n,t,y)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400  enddo
 500  return
      end

      subroutine kleap(kappa,t,f,y,v,d,n,m,ntreeh,k,iflag)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      real gengam
      dimension y(n),f(m,ntreeh),p(m),v(n,m),d(n,m)
      dimension k(m),ichangey(n)
      external gengam,genmul,fprob
c=========================================
c Determine time interval and update time
c=========================================
      fsum=f(1,ntreeh)
      if(fsum.gt.0.0d0)then
         tstep=gengam(real(fsum),real(kappa))
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
      call genmul(kappa,real(p),m,k)
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
               f(j,1)=fprob(j,n,t,y)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400  enddo
 500  return
      end
	  
      subroutine tauleap(tau,t,f,y,v,d,n,m,ntreeh,k,iflag)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension y(n),f(m,ntreeh),k(m),v(n,m),d(n,m),ichangey(n)
      external poidev,fprob
c======================================================================
c Generate Poisson random variables (number of event firings in t+tau)
c======================================================================
      fsum=f(1,ntreeh)
      if(fsum.gt.0.0d0)then
         do j=1,m
            k(j)=ignpoi(real(f(j,1)*tau))
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
               f(j,1)=fprob(j,n,t,y)
               diff=f(j,1)-fold
               jdum=int((j+1)/2)
               do itree=2,ntreeh
                  f(jdum,itree)=f(jdum,itree)+diff
                  jdum=int((jdum+1)/2)
               enddo
               goto 400
            endif
         enddo
 400  enddo
 500  return
      end
