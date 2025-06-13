c
c returns Born (nlo=0) or 
c finite part of virtual contributions (NOT sum of Born + virtuals!) 
c
c---------------------------------------------------------------------------
c
      subroutine qq_ee_vv_custom(pbar,sign, nlo, L,k,id_beam,ans)
      implicit none
c
c	Last modified for POWHEG by Barbara Jaeger: 2024 Jan.
c
C  this routine calculates the matrix elements**2 for quark lepton scattering
c  (gamma or Z or W exchange) in the t-channel:
C
C        l1 q3    ---->   l2 q4 (l can be charged lepton or neutrino)
c
C  Crossing related processes can be computed as well. 
c  Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production is NOT implemented.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return qqll = 2Re(M_Born^* M_virt)
c          NLO = 0       return qqll = |M_born|^2   etc.
c
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)
      include 'pwhg_st.h'

      include 'global.inc'
      include "coupl_custom.inc"
      
      include 'vtype.h'
      include 'pwhg_flg.h'
c      
      integer mu 

c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c alfas, scales etc
      include 'scales.inc'
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+nv), musq
      double precision ans
      double precision res(3),resv(3)
      double precision  p(0:3,4+nv), p21(0:4), p43(0:4)
      integer  sign(4+nv), nlo, i, j, jj, k, kk, id,
     1     isig, isig1, isig3,h,kl
      integer id_beam

      integer  ifl(4,3), js1, js3, L, is1, is3
      integer ipol
      double complex mat(3,-1:1,-1:1,3) 
      double complex matv(3,-1:1,-1:1,3)
      double complex mm(3,-1:1,-1:1), 
     1               mv12(3,-1:1,-1:1),
     1               mv34(3,-1:1,-1:1)
      double complex mvv,maa,mzz, mww
      double complex den_a,den_z,den_w
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) 
      double complex fac
      double complex dotcc 
      external dotcc 

      logical nc_type,cc_type
      logical linit
      data linit /.true./ 
      
      save ifl,  linit 
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ii,ll
c
c variables for powheg:
      double precision q2_up,q2_lo,rup,rlo,lrup,lrlo
      double precision cvirtl   

c identify "bad" points (low qsq):
      logical qdamp   

c select specific contributions:
c      logical zonly,aonly
c      parameter (zonly=.false.,aonly=.false.)
c
c variables for virtual corrections   
c
      double precision c2,c2o4pi     !,pi2o3, cvirtc
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo

      lnlo = NLO.ne.0    ! include some virtual stuff if T

c  ---------------------------------------------------------------------
c ...IB: polarization variable
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise  
      ipol = merge(1, 0, flg_polarised_dis)
c  ---------------------------------------------------------------------
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
c	
      if (k.lt.1) stop 'this value of k is not allowed: k = '
      if (k.gt.3) stop 'this value of k is not allowed: k = '

c fix strong coupling gs**2 for the 2 quark lines:
      als = st_alpha
c
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
c
         print*,'Born/virtual amplitudes:'

         if (id_beam.eq.1) then ! charged lepton beam
            kl = 1              ! ee-uu
            ifl(1,kl) = 2
            ifl(2,kl) = 2
            ifl(3,kl) = 3       
            ifl(4,kl) = 3
            kl = 2              ! ee-dd
            ifl(1,kl) = 2
            ifl(2,kl) = 2
            ifl(3,kl) = 4
            ifl(4,kl) = 4
            kl = 3              ! ev-ud
            ifl(1,kl) = 2
            ifl(2,kl) = 1
            ifl(3,kl) = 3
            ifl(4,kl) = 4       

         else                   ! neutrino beam
            kl = 1              ! vv-uu
            ifl(1,kl) = 1
            ifl(2,kl) = 1
            ifl(3,kl) = 3       !fermion=lepton type
            ifl(4,kl) = 3
            kl = 2              ! vv-dd
            ifl(1,kl) = 1
            ifl(2,kl) = 1
            ifl(3,kl) = 4
            ifl(4,kl) = 4
            kl = 3              ! ve-ud
            ifl(1,kl) = 1
            ifl(2,kl) = 2
            ifl(3,kl) = 3
            ifl(4,kl) = 4       
         endif !id_beam
         
      endif !linit

      if (k.le.2) then
         nc_type = .true.
         cc_type = .false.
      else
         cc_type = .true.
         nc_type = .false.
      endif

      mat = 0d0
      matv = 0d0
c      do kl = 1,3
c         do isig1 = -1,1,2
c            do isig3 = -1,1,2
c               do i = 1,3 (gamma,Z,W)
c                  mat(kl,isig1,isig3,i) = 0
c                  matv(k,isig1,isig3,i) = 0
c               enddo !i
c            enddo 
c         enddo !isig
c      enddo
c
c identify fermion line sign factors
c
      is1 = sign(1)
      is3 = sign(3)
      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1
c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 
         p21(mu) =   p(mu,2) - p(mu,1)
         p43(mu) =   p(mu,4) - p(mu,3)	 
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2

      qdamp = .false.

      if (abs(p21(4)).lt.qsqAmin.or.abs(p43(4)).lt.qsqAmin) 
     &     qdamp = .true.

c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2)) 
c
c -------------------------------------------------------------------
c
c contract currents to get t-channel diagrams
c
c use Diract equation: ubar(p2)(p2-p1)_slash u(p1) = 0 etc. 
c to disregard q.mu q.nu term in V propagator
c
c
      if (nc_type) then

         den_z = dcmplx(p43(4)-xm2(2),xmg(2))
         den_a = dcmplx(p43(4),0d0)

      do isig1 = -1,1,2
         do isig3 = -1,1,2
c
          mvv = dotcc(jqq(0,isig1,1),jqq(0,isig3,2))

          maa = mvv/den_a
          mzz = mvv/den_z

c          ! photon exchange:
             mat(k,isig1,isig3,1) = 
     4             maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
c          ! Z exchange:    
             mat(k,isig1,isig3,2) = 
     4             mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)

             matv(k,isig1,isig3,1) = (0d0,0d0)
             matv(k,isig1,isig3,2) = (0d0,0d0)
	    
         enddo !isig3
      enddo !isig1
      
c----------------------      
    
      else ! cc_type  

         den_w = dcmplx(p43(4)-xm2(3),xmg(3))

c W exchange:
         if (k.eq.3) then 
            mww = dotcc(jqq(0,-1,1),jqq(0,-1,2))/den_w
            mat(k,-1,-1,3) = mww*clr(ifl(1,k),3,-1)*clr(ifl(3,k),3,-1)
            matv(k,-1,-1,3) = (0d0,0d0)
         endif

      endif !nc/cc
c
c
c -------------------------------------------------------------------------
c -------------------------------------------------------------------------

999     continue

        if(zonly) mat(:,:,:,1) = 0d0 ! no gamma exchange
        if(aonly) mat(:,:,:,2) = 0d0 ! no Z exchange
     
c sum the graphs, square them and map them onto uucc, uuss etc.

c      do k = 1,3
         res(k) = 0
         resv(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,3
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                  mat(k,isig1,isig3,i)
               enddo !i
c       
c               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
c     &                         + dimag(mm(k,isig1,isig3))**2


c ...IB: polarized cross section: 
c ...    (note that for antiparticles isig = -physical helicity      
c ...    --> added extra factor is1*is3)
c
            if(ipol.eq.0) then

               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2

            elseif(ipol.eq.1) then                

               res(k) = res(k) + is1*is3*( 
     &                         isig1*isig3*dreal(mm(k,isig1,isig3))**2
     &                        +isig1*isig3*dimag(mm(k,isig1,isig3))**2)

            else 

               write(6,*) 'ipol can be only 0 or 1, aborting'
               STOP

            endif


		mv12(k,isig1,isig3) = 0d0
		mv34(k,isig1,isig3) = 0d0

                if (lnlo.and.(.not.qdamp)) then
                  
c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the factor pi^2/3+9/2 for the born term is 
c  after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
c
c  comply with POWHEG normalization:
c  vertex corrections in dreg:
c
c   V^nu = (4*Pi)^ep * Gamma(1+ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7)*B^nu
c
c         = (4*Pi)^ep / Gamma(1-ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7-Pi^2/3)*B^nu
c
c     The factor  (4*Pi)^ep/Gamma(1-ep) IS NOT RETURNED by this subroutine
c     and it's thought as factorized in front of the real counterterms too.
c
c     if (4*Pi)^ep / Gamma(1-ep) is collected in front then cvirt:
c     set parameter in dred (cvirtl = -8d0)
c
c     squared momentum of the weak boson connected with the upper line
          q2_up = p21(4)
c     squared momentum of the weak boson connected with the lower line
          q2_lo = p43(4)

c          rup = st_muren2/(-q2_up)
c          if (rup.lt.0d0.and..not.qdamp) then
c              write(*,*) 'Error in virtuals: q2_up should be < 0!!'
c              write(*,*) 'q2_up = ', q2_up
c              write(*,*) 'st_muren2 = ', st_muren2
c              ! set to dummy value:
c              	rup = 1d0
c          endif      
          rlo = st_muren2/(-q2_lo)
          if (rlo.lt.0d0.and..not.qdamp) then
              write(*,*) 'Error in setvirtual: q2_lo should be < 0!!'
              write(*,*) 'q2_lo = ', q2_lo
              write(*,*) 'st_muren2 = ', st_muren2
              ! set to dummy value:
              rlo = 1d0
              stop
          endif      
c          lrup = log(rup)
          lrlo = log(rlo) 

          cvirtl = -8d0
c
                  if (nlo.gt.0) then
c                     mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
c     1                ( mv12(k,isig1,isig3) + 
c     2                  mm(k,isig1,isig3)*(-lrup**2-3*lrup+cvirtl) )
c
                     if(ipol.eq.0) then
                        mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                  ( mv34(k,isig1,isig3) + 
     2                  mm(k,isig1,isig3)*(-lrlo**2-3*lrlo+cvirtl) )
c
c ...IB: polarized cross section:
c ...    (note that for antiparticles isig = -physical helicity      
c ...    --> added extra factor is1*is3)
c     
                     elseif(ipol.eq.1) then
                        mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                  ( mv34(k,isig1,isig3) + 
     2                  (is1*is3)*isig1*isig3*
     3                  (mm(k,isig1,isig3)*(-lrlo**2-3*lrlo+cvirtl) ))
                     else 
                        write(6,*) 'set polarization in powheg.input'
                        stop
                     endif
                  else
c                     mv12(k,isig1,isig3) = 
c     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                     mv34(k,isig1,isig3) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3) * conjg( mv34(k,isig1,isig3) )  )

c                  resv(k) = resv(k) + 2*dreal(
c     1                 mm(k,isig1,isig3) * conjg( mv12(k,isig1,isig3) )  )

               endif !lnlo
           enddo !isig3
         enddo !isig1 

         if (nlo.eq.0) then
            res(k) = res(k)*3d0
         elseif (nlo.gt.0) then
            if (qdamp) then
               !fakevirt:
                res(k) = res(k)*3d0   ! Born-type
                res(k) = res(k)*0.2d0 ! 'fakevirt' factor
            else !no damping  
c virt only (without Born):
               res(k) = (resv(k))*3d0 ! 3 is the color sum factor
            endif !qdamp
         else !nlo.lt.0
            if (qdamp) then
               !fakevirt:
                res(k) = res(k)*3d0   ! Born-type
                res(k) = res(k)*0.02d0 ! 'fakevirt' for boxes, pentagons
            else !no damping  
c virt only (without Born):
               res(k) = resv(k)*3d0 ! 3 is the color sum factor
            endif !qdamp
         endif !lnlo
c      enddo !k

           
cc eliminate processes with photon virtuality below cutoff
      if (qdamp) then 
          res(k) = res(k)*1d-20
      endif

      ans = res(k)
      return
      end


