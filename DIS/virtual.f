      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegbornexternal)
      real * 8 p(0:3,nlegs)
      integer vflav(nlegs)
      real * 8 born,virtual
      integer j,k,mu,nu

      real *8 powheginput
      external powheginput 

      logical, save :: firsttime = .true. 
      logical, save :: secondtime = .false. 

      integer fakevirt
      save fakevirt 
      integer, save :: flst_counter = 0
      integer nmomset
      parameter (nmomset=10)

      
      real *8 q2cutB, q2cutR, q2cutA
      common /q2cut/q2cutB, q2cutR, q2cutA

      real * 8, external :: dotp
cccccccccccc

      if(2d0*dotp(p(:,1), p(:,3))<q2cutB) return

      if (kn_jacborn.eq.0d0) then
         virtual = 1d-30
         return
      endif

      if (firsttime) then
         fakevirt=powheginput("#fakevirt")
         if (fakevirt == 1) then 
            write(*,*) 'WARNING: Using fakevirt !'
            firsttime=.false.
         endif
         if(powheginput("#parallelstage")<2) then
            fakevirt = 1        ! First point ALWAYS run with fakevirt. This ensures that smartsig works properly
                                ! Only doing it for the first stage, otherwise problems during event generation
            flst_counter=flst_counter+1
            if(flst_counter.ge.(flst_nborn*nmomset)) then
               firsttime = .false.
               secondtime=.true.
            endif
         else
            firsttime = .false.
         endif
      elseif (secondtime) then 
         fakevirt=powheginput("#fakevirt")
         secondtime = .false.
      endif

      
      if(fakevirt.eq.1) then  

         call compborn_eq(p,vflav,born) 
         virtual = 0.2d0*born
      
      else

c numbering of momenta is e(1)q(2)->e(3)q(4)
c
         call compvirt_eq(p,vflav,virtual) 
c         
c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
         virtual = virtual/(st_alpha/(2d0*pi))
         
c~          call compborn_eq(p,vflav,born)
c~          virtual = -8d0*4d0/3d0*born
         
      endif
      
      return      
      end
c     
c==================================================
c
      subroutine compvirt_eq(pin,bflav,virtual)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
c
      integer nlegs,nf
      parameter (nlegs=nlegbornexternal)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 virtual
c
c vbfnlo stuff:
      include 'global.inc'
      integer nlo
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,4+nv), polcol
      real*8 res

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,j,mu,nu
      integer FSIGN(4+nv),physToDiag(4)
      
      integer ftype(nlegs)
      integer k,id_beam

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nlo = 1 ! NLO
c one quark, one electron in initial state:
      polcol = 1d0/(4d0*N)

      ftype = 1
      
      virtual = 0d0

      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         
         p(mu,3) = pin(mu,3)
         p(mu,4) = pin(mu,4) 

         p(mu,5) = 0d0   

      enddo ! mu
      ! Particle ordering of vbfnlo is (powheg is 1234)
      ! e(1) q/qb(3) -> e(2) q/qb(4)


      fsign(1) = sign(1,bflav(1))
      fsign(2) = sign(1,bflav(3))
      fsign(3) = sign(1,bflav(2))
      fsign(4) = sign(1,bflav(4))

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
      
      if (bflav(2).gt.0) then       
            physToDiag(2)=3
            physToDiag(4)=4
      else
            physToDiag(2)=4
            physToDiag(4)=3
      endif
      
C*****************  end of process evaluation  **********************

c get the amplitude squared:
      do mu = 0,3
         do i = 1,4
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo	 

c
c eq->eq or eq->vq (for v-beam: vq->vq or vq->eq):
c
c NC:
c eueu-type(channel_type=1)
c eded-type(channel_type=2)
c
c CC:
c euvd-type(channel_type=3)

c use value of bflav to select appropriate uucc, uuss etc.
      if(mod(abs(bflav(2)),2).eq.0) ftype(2) = 2 ! fermion1 = up-type 
      if(mod(abs(bflav(4)),2).eq.0) ftype(4) = 2 ! fermion3 = up-type 

      if ((ftype(2).eq.ftype(4))) then !nc
         k = -ftype(2)+3
      else !cc 
         k = 3
      endif !nc/cc  

      id_beam=1                 !lepton beam
      if ((abs(bflav(1)).eq.12).or.(abs(bflav(1)).eq.14).or.
     &    (abs(bflav(1)).eq.16)) id_beam=0 !v beam
      
      call qq_ee_vv_custom(pbar,fsign,1,1,k,id_beam,res)

c     Sum over polarization for neutrino-induced is 1/2 (only left handed nu),
c     for charged leptons it is 1/4:                                      
      polcol = polcol*dble(2-id_beam)
      
      virtual = res*polcol

      return
      end




      
!     Silvia's simple routines below without the Z
!     subroutine setvirtual(p,bflav,virtual)
!      implicit none
!      include 'pwhg_math.h'
!      include 'nlegborn.h'
!      include 'pwhg_flst.h'
!      include 'pwhg_st.h'
!      include 'pwhg_kn.h'
!
!      integer nlegs
!      parameter (nlegs=nlegbornexternal)
!      real * 8 p(0:3,nlegs)
!      integer bflav(nlegs)
!      real * 8 virtual, q2, L
!      real * 8, external ::dotp
!
!      q2 = 2d0 *dotp(p(:,1), p(:,3))
!      L = log(st_muren2/q2)
!
!      call setborn_onlyamp(p,bflav,virtual)
!      virtual = virtual * (-8d0 -3d0 * L - L**2) * cf 
!      end
