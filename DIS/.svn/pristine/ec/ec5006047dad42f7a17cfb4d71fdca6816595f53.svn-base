      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegbornexternal)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs), myborn, s, t, u, r
      real * 8 myp(0:3,nlegs)
      integer bflav(nlegs), mybflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf,born_old
      integer j,k,mu,nu
      real *8, external :: dotp

      real *8 q2cutB, q2cutR, q2cutA
      common /q2cut/q2cutB, q2cutR, q2cutA

      
c Colour factors for colour-correlated Born amplitudes;
c Rule from 2.98 in FNO2007, leads to B_i j=B*(C_i+C_j-C_k)/2,
c where k#i,j

c compute Born amplitude squared for e q -> e q: 
c     (assume particles 1 and 3 are leptons, 2 and 4 are quarks)

      call compborn_eq(p,bflav,born)

      bmunu = 0d0      
      do j=1,nlegs
         if(abs(bflav(j)).le.6) then
            do k=j+1,nlegs
               if (((j.eq.2).and.(k.eq.4))) then   
                  colcf = cf            
               else
                  colcf = 0
               endif
               bornjk(j,k)=born*colcf
               bornjk(k,j)=bornjk(j,k)
            enddo
         endif
      enddo

      end
c     
c==================================================
c
      subroutine compborn_eq(pin,bflav,born)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
c
      integer nlegs,nf
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 born
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
      nlo = 0 ! LO
c      polcol =1d0/(4d0*N**2)
c one quark, one electron in initial state:
      polcol = 1d0/(4d0*N)

      ftype = 1
      
      born = 0d0

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
      endif                     !nc/cc

      id_beam=1 !lepton beam
      if ((abs(bflav(1)).eq.12).or.(abs(bflav(1)).eq.14).or.
     &    (abs(bflav(1)).eq.16)) id_beam=0 !v beam
      
      call qq_ee_vv_custom(pbar,fsign,0,1,k,id_beam,res)

c     Sum over polarization for neutrino-induced is 1/2 (only left handed nu),
c     for charged leptons it is 1/4:                                      
      polcol = polcol*dble(2-id_beam)
      
      born = res*polcol

      return
      end


      
!     Silvia's simple amplitudes without Z below
!     subroutine setborn(p,bflav,born,bornjk,bmunu)
!      implicit none
!      include 'pwhg_math.h'
!      include 'nlegborn.h'
!      include 'pwhg_flst.h'
!      include 'PhysPars.h'
!      integer nlegs
!      parameter (nlegs=nlegbornexternal)
!      integer bflav(nlegs)
!      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs), bmunu(0:3, nlegs, nlegs), born
!      real *8, external :: dotp
!
!      call  setborn_onlyamp(p,bflav,born)
!      
!      bmunu  = 0d0
!      
!      bornjk = 0d0
!      bornjk(2,4) = born * cf
!      bornjk(4,2) = born * cf
!      end
!
      subroutine setborn_onlyamp(p,bflav,born)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegbornexternal)
      real * 8 p(0:3,nlegs), born, s, t, u
      integer bflav(nlegs)
      real *8, external :: dotp  


      s =2d0*dotp(p(:,1),p(:,2))
      t=-2d0*dotp(p(:,1),p(:,3))
      u=-2d0*dotp(p(:,2),p(:,3))

      born = ph_unit_e**4 * (u**2+s**2)/t**2 *2d0
      if(mod(bflav(2),2) ==0) then
         born = born*(2d0/3d0)**2
      else
         born = born /(3d0)**2
      endif
      end

      


c==================================================
c
      subroutine borncolour_lh
      implicit none
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structures, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface

      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer i

C     -- neutral particles
      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,3)=0
      icolup(2,3)=0

c     -- colored particles
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,4)=0
      icolup(2,4)=0

      if(idup(2).gt.0) then
         icolup(1,2)=501
         icolup(2,2)=0
      else
         icolup(1,2)=0
         icolup(2,2)=501
      endif
      
      do i=1,2
         icolup(i,4)=icolup(i,2)
      enddo
      end
c

cccccccccccccc

      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
      
      call lhefinitemasses
      end
