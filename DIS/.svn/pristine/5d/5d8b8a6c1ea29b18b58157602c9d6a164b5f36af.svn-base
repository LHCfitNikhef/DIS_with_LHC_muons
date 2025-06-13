c
      subroutine setreal(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'

      integer nleg
      parameter (nleg=nlegbornexternal+1)
      real * 8 p(0:3,nleg)
      integer fermion_flav(nleg)
      real * 8 amp2

      real * 8, external ::dotp
      real *8 q2cutB, q2cutR, q2cutA
      common /q2cut/q2cutB, q2cutR, q2cutA
	 
c----------------------------------------------------
c
      amp2 = 0d0
      if(2d0*dotp(p(:,1),p(:,3))< q2cutR) return

      
      if (kn_jacborn.eq.0d0) then 
         amp2 = 1d-30
         return
      endif
      
      call compreal_eq(p,fermion_flav,amp2)

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compreal_eq(pin,bflav,amp2)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
c
      integer nlegs
      parameter (nlegs=nlegbornexternal + 1)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 amp2 
c
c vbfnlo stuff:
      include 'global.inc'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,5+nv),qbar(0:4)
      real*8 polcol,polcolq,polcolg
      real*8 res(3)

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu
      integer FSIGN(4+nv),gsign,physToDiag(5)
      
      real*8 nans(2,2,3),cans(2,3)
      logical nc_type
      integer k,icc,id_beam

      integer ftype(nlegs)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      polcol = 0d0
c      polcolq = 1d0/(4d0*N**2)
c      polcolg = 1d0/(4d0*N*(N**2-1))
c one quark, one electron in initial state:
      polcolq = 1d0/(4d0*N)
      polcolg = 1d0/(4d0*(N**2-1))

      ftype = 0
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   

         if (bflav(2).ne.0) then ! final-state gluon 
            if (bflav(5).eq.0) then
               p(mu,3) = pin(mu,3)
               p(mu,4) = pin(mu,4) 
               p(mu,5) = pin(mu,5)     
            elseif (bflav(4).eq.0) then
               p(mu,3) = pin(mu,3)
               p(mu,4) = pin(mu,5) 
               p(mu,5) = pin(mu,4)  
            endif
         else   ! initial-state gluon 
            p(mu,3) = pin(mu,3)
            p(mu,4) = pin(mu,4) 
            p(mu,5) = pin(mu,5) 
         endif ! fin/in state gluon
      enddo ! mu

!     The e-/e+ is always in position 1 (incoming) and 2 (outgoing) in vbfnlo 
!     and always in position 1 (incoming) and 2 (outgoing) in powheg. If we use 
!     anti-electrons then wecross them (ie a n incoming positiron becomes an 
!     outgoing electron). Same for quarks.
      if(bflav(1).gt.0) then
            physToDiag(1)=1             
            physToDiag(3)=2                  
      else
            physToDiag(1)=2             
            physToDiag(3)=1
      endif

      fsign(1) = sign(1,bflav(1))
      fsign(2) = sign(1,bflav(1))   

      if (bflav(2).gt.0.and.bflav(4).gt.0) then

C*******************  e1 q3 ---> e2 q4 g   **********************

      physToDiag(2)=3
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol = polcolq

c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(2)),2)         
      ftype(4) = 2-mod(abs(bflav(4)),2)

      elseif (bflav(2).lt.0.and.bflav(4).lt.0) then
c            
C******************* e1 qbar4 ---> e2 qbar3 g   **********************
      
      physToDiag(2)=4
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq
      
c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(2)),2)         
      ftype(4) = 2-mod(abs(bflav(4)),2)      

      elseif (bflav(2).gt.0.and.bflav(5).gt.0) then

C*******************  e1 q3 ---> e2 g q  **********************

      physToDiag(2)=3
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol = polcolq

c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(2)),2)         
      ftype(4) = 2-mod(abs(bflav(5)),2)

      elseif (bflav(2).lt.0.and.bflav(4).eq.0) then
c            
C******************* e1 qbar4 ---> e2 g qbar3   **********************
      
      physToDiag(2)=4
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq
      
c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(2)),2)         
      ftype(4) = 2-mod(abs(bflav(5)),2)      

      elseif (bflav(2).eq.0.and.bflav(4).gt.0) then
      
C*******************  e g ---> e q qb   **********************
      
      physToDiag(2)=5
      physToDiag(4)=4
      physToDiag(5)=3 
c
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg
c
c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(5)),2)         
      ftype(4) = 2-mod(abs(bflav(4)),2)   

      elseif (bflav(2).eq.0.and. bflav(4).lt.0) then
      
C*******************  e g ---> e qb q   **********************

      physToDiag(2)=5            
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ftype(2) = 2-mod(abs(bflav(4)),2)         
      ftype(4) = 2-mod(abs(bflav(5)),2)   

      else
         
         print*,'this real-flavor combination is not implemented'
         print*,'rflav=',bflav
         stop

      endif

c**************************************

      
!     Determine if it is CC or NC based on flavour structure of quarks
      if ((ftype(2).eq.ftype(4))) then
         k = -ftype(2)+3
      else !cc 
         k = 3
      endif                   

      id_beam=1 !lepton beam
      if ((abs(bflav(1)).eq.12).or.(abs(bflav(1)).eq.14).or.
     &    (abs(bflav(1)).eq.16)) id_beam=0 !v beam

      
C*****************  end of process evaluation  **********************

      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo 
      qbar(4) = 0d0

         
c k-ordering is: eueu,eded,euvd
      call qqj_ee_vv_custom(pbar,fsign,qbar,gsign,k,id_beam,res(1))

c     Sum over polarization for neutrino-induced is 1/2 (only left handed nu),
c     for charged leptons it is 1/4:                                      
      polcol = polcol*dble(2-id_beam)
            
      amp2 = res(1)*polcol

      return
      end


      
c     Silvia's simple routines without Z below
!      subroutine setreal(pin,rflav,amp2)
!      implicit none
!      include 'nlegborn.h'
!      include 'pwhg_math.h'
!      include 'pwhg_st.h'
!      include 'pwhg_kn.h'
!      include 'PhysPars.h'
!
!      integer nleg
!      parameter (nleg=nlegbornexternal+1)
!      real * 8 pin(0:3,nleg)
!      integer rflav(nleg)
!      real * 8 amp2, myamp2
!
!      real * 8 p12, p15, p24, p25, p45, q2
!      real * 8, external :: dotp
!      real * 8 p(0:3, nleg)
!
!      p=pin
!      ! Swap leg 2 and 5 if there is a gluon in the initial state
!      if(rflav(2) .eq. 0) then
!         p(:, 5) = -pin(:,2)
!         p(:, 2) = -pin(:,5)
!      endif
!      
!
!      p12 = dotp(p(:,1),p(:,2))
!      p15 = dotp(p(:,1),p(:,5))
!      p24 = dotp(p(:,2),p(:,4))
!      p25 = dotp(p(:,2),p(:,5))
!      p45 = dotp(p(:,4),p(:,5))
!
!      q2 = 2d0* dotp(p(:,1),p(:,3))
!
!      ! Amplitude for e q (qbar) -> e q (qbar) g
!      if(p45.ne.0d0) then
!         myamp2=        (2*(4*p12**2 + 2*p15**2 + 2*p24**2 + 4*p24*p25 + 3
!     $        *p25**2 + 2*p15*(p24 + 2*p25 - p45) - 2*p12*(2*p15 + 2*p24 +
!     $        3*p25 - p45) - 2*p24*p45 - 2*p25*p45 + p45**2))/(p25*(p24 +
!     $        p25 - p45)*p45)
!      else
!         myamp2 = 0d0
!      endif
!      myamp2 = myamp2 * CF * ph_unit_e**4 * 8 *pi**2
!
!
!      if(rflav(4).eq. 0) then
!         print *, "Stopping as the code assumes 4 to be a fermion"
!         print *, rflav
!         call pwhg_exit(-1)
!      endif
!      
!      ! Determine the charge of the quark for the QED coupling
!      if(mod(rflav(4),2) ==0) then
!         myamp2 = myamp2*(2d0/3d0)**2
!      else
!         myamp2 = myamp2 /(3d0)**2
!      endif
!
!
!!     Incoming gluon channel
!!     Minus sign from pertmuting a fermin line
!!     Factor Tf/CF as we need to divide by the number of colours for the incoming 
!      if(rflav(2)  == 0) then
!         myamp2 = -myamp2/CF*TF
!      endif
!
!      amp2=myamp2
!      
!      end


c There are no regulars
      subroutine regularcolour_lh
      implicit none
      end
