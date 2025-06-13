
      subroutine qqj_ee_vv_custom(pbar,sign,qbar,gsign, k,id_beam,ans)
      implicit none
c
c	Last modified for POWHEG:  Oct. 2019
C
C this routine calculates the matrix elements**2 for quark electron scattering
C
C        e1 q3    ---->   e2 q4 g 
c
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production is NOT implemented.
c
C  This code includes only real emission contributions, i.e.
c
c      return uuee = |M_real|^2   etc.
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper (12) line (j=2) or lower (34) line (j=3)
c	l is the gluon polarization in the kartesian basis (l=1,2)
c         l=0 stands for building blocks without gluon emission
c	k is the process ID (1:uuee,2:ddee,3:udev)
c	isig1/isig3 are the helicities of parton1/3 
c
c---------------------------------------------------------------------
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
      double complex den_w,den_z,den_a
c
      double precision q12(0:4,3),q34(0:4,3)
      
c electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c
c alfas, scales etc
      include 'scales.inc'
c
c variables for the main part of the program
c
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,4+nv),qbar(0:4),musq
      double precision  res(6,2:3)
      double precision ans(3)     
      integer id_beam
      
      integer kl,h,mu
      double precision  p(0:3,4+nv),q(0:4), 
     1                  pq(0:4,4)

      integer  sign(4+nv),gsign, i, j,jj, k, id,
     1         isig, isig1, isig3
      integer  ifl(4,6), js1, js3, is1, is3
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      integer ipol
      logical jlog2,jlog3
      double complex mat(6,-1:1,-1:1,2:3,0:2,3)
      double complex mat_ghe(6,-1:1,-1:1,2:3,0:2,3)      
      double complex mm(6,-1:1,-1:1,2:3,2)
      double complex mvv,maa,mzz,mww
      double precision eps(0:3,2) ! g in kartesian basis 
c      double complex eps_temp(0:3,2)
      double complex psi(2,-1:1,4),jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex dotcc,dotrc
      external dotcc,dotrc
      integer lh

      double complex zm2i(2:3),ezz1,ezz2,ezz
c
      logical nc_type,cc_type
c
      logical linit
      save linit,ifl
      data linit /.true./

c select specific contributions:
c      logical zonly,aonly
c      parameter (zonly=.false.,aonly=.true.)
c      parameter (aonly=.false.,zonly=.true.)
c       parameter (zonly=.false.,aonly=.false.)

      logical qdamp
      parameter (qdamp=.false.)
c
c  ---------------------------------------------------------------------
c ...IB: polarization variable
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise  
      ipol = merge(1, 0, flg_polarised_dis)
c
c  ---------------------------------------------------------------------
c
c initialize & precompute stuff needed below:
c
c fix strong coupling gs**2 for the 2 quark lines:
c
      als = st_alpha

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
	  
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
         print*,'real-emission amplitudes:'

         if (zonly) print*,'NC: only Z exchange diagrams (no photon)'
         if (aonly) print*,'NC: only photon exchange diagrams (no Z)'
         if (qdamp) then 
            print*
            print*, 'minimum virtuality for t-channel boson exchange'
            print*, 'qsqA_min = ',qsqAmin,'GeV**2'
            print*,'damping factor of 1d-20 below qsqAmin '
         endif !qdamp
c
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
            ifl(3,kl) = 3       
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

      if (gsign.eq.1) then	!final state gluon
        jlog2 = .true.		! can couple to upper/lower line
	jlog3 = .true. 
      	jmin = 2
	jmax = 3
      else		       !initial state gluon -> only:
	 if (sign(1).ne.sign(2)) then  !gluon from upper line
	       jlog2 = .true.
	       jlog3 = .false. 
	       jmin = 2
	       jmax = 2
	 else			       !gluon from lower line
	       jlog2 = .false.
	       jlog3 = .true. 
	       jmin = 3
	       jmax = 3
       endif
      endif

      mat = 0d0
c      do kl = 1,6
c         do isig1 = -1,1,2
c            do isig3 = -1,1,2
c	       do j = 2,3
c	          do l = 0,2
c                  	mat(kl,isig1,isig3,j,l,i=1:3)  = 0d0
c		  enddo	!l
c               enddo !j
c            enddo
c         enddo !isig
c      enddo
c
c identify fermion line sign factors (for 1 3 -> 2 4 etc.)
c
c      is1 = sign(1)
c      is3 = sign(3)
c      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
c      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1

c fix is1 such that is1 = +1 for   q1 ->  q2 g 
c		    is1 = -1 for  q1b -> q2b g 
c		    is1 =  0 for    g -> q1b q2
c
c      for is3:     is3 = +1 for   q3 ->  q4 g 
c		    is3 = -1 for  q3b -> q4b g 
c		    is3 =  0 for    g -> q3b q4

	is1 = (sign(1)+sign(2))/2  
	is3 = (sign(3)+sign(4))/2
		
c (is1,is3 are fixed here and don't change throughout this run of the program)
c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 q(mu) = qbar(mu)*gsign	 
      enddo
      q(4)   = 0d0

c  ---------------------------------------------------------------------
c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)

c get the f-fbar currents (with no gluon radiation) 
c	J21^mu=jqq(mu,isig1,1,is1,0), J43^mu=jqq(mu,isig3,2,is3,0) 
c
        call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),
     #						jqq(0,-1,1,is1,0))      
        call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),
     #						jqq(0,-1,2,is3,0))

c  Get the gluon polarization vector
      do l = 1,2	! 2 gluon polarizations
         call polvec(qbar,l,eps(0,l))
c         call helvec(qbar,gsign,l,eps_tmp(0,-2*l+3)) !IB debug        
	 
c for gauge check:
c        if (lgauge) then ! contract amplitude with q rather than eps(q)
c	 do mu = 0,3
c	 	eps(mu,l) = qbar(mu)
c	 enddo		 
c	 endif
	 	 
         do isig = -1,1,2	! fermion helicity 
 
c     NOTE for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise: 

c upper line emission:            
c            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
c     $           braketg(1,isig,1,l),pq(0,1))      
c            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
c     $           braketg(1,isig,2,l),pq(0,2))      
     
c lower line emissions:
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,l),
     $           braketg(1,isig,3,l),pq(0,3))      
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,l),
     $           braketg(1,isig,4,l),pq(0,4))    

         enddo !isig
      enddo !l
       
c     Get the f-fbar currents with one gluon radiated from the
c     current line:
c
c	gluon from upper line:
      do l = 1, 2 ! gluon polarizations
c         call curr6(1,psi(1,-1,2),p(0,2),braketg(1,-1,1,l),pq(0,1),jh1)	
c         call curr6(1,braketg(1,-1,2,l),pq(0,2),psi(1,-1,1),p(0,1),jh2)	
c         do isig = -1,1,2 ! fermion helicity
c            do mu = 0,5
c	       jqq(mu,isig,1,is1,l) = jh1(mu,isig) + jh2(mu,isig)
c            enddo
c         enddo
c never have gluon from upper line:
         jqq(:,:,1,is1,l) = 0d0 
         
c gluon from lower line:
c         jqq(:,:,2,is3,l) = 0d0 
         call curr6(1,psi(1,-1,4),p(0,4),braketg(1,-1,3,l),pq(0,3),jh1)	
         call curr6(1,braketg(1,-1,4,l),pq(0,4),psi(1,-1,3),p(0,3),jh2)	
         do isig = -1,1,2
            do mu = 0,5
               jqq(mu,isig,2,is3,l) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo

      enddo !l
      

      do mu = 0,3
c         q12(mu,1) = p(mu,1)-q(mu)-p(mu,2)
         q12(mu,1) = p(mu,1)-p(mu,2)
         q34(mu,1) = p(mu,3)-q(mu)-p(mu,4)         
c adapted for POWHEG:
            q12(mu,2) = q12(mu,1)
            q12(mu,3) = p(mu,1)-p(mu,2)
            q34(mu,2) = p(mu,3)-p(mu,4)
            q34(mu,3) = q34(mu,1)
      enddo
      do L = 1,3
         q12(4,L) = q12(0,L)**2-q12(1,L)**2-q12(2,L)**2-q12(3,L)**2
         q34(4,L) = q34(0,L)**2-q34(1,L)**2-q34(2,L)**2-q34(3,L)**2
      enddo
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c contract currents to get t-channel diagrams
c
c use Diract equation: ubar(p2)(p2-p1)_slash u(p1) = 0 etc. 
c to disregard q.mu q.nu term in V propagator
c
c neutral current first:
c
      if (nc_type) then    

         den_z = dcmplx(q12(4,3)-xm2(2),xmg(2))
         den_a = dcmplx(q12(4,3),0d0)

c         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))

      do l = 1,2	! gluon polarization
        do isig1 = -1,1,2  ! fermion1 helicity
          do isig3 = -1,1,2  ! fermion3 helicity
	 
           if (jlog3) then
	   j = 3 ! g from lower line
	   
           mvv = dotcc(jqq(0,isig1,1,is1,0),jqq(0,isig3,2,is3,l))
           
c extra gauge term for Z: (q12.J2) (q12.J1)/msq/(prop factor)
c not needed, because vanishes via Dirac equation (lepton line):
c           ezz1 = dotrc(q12(0,3),jqq(0,isig1,1,is1,l))
c           ezz2 = dotrc(q12(0,3),jqq(0,isig3,2,is3,0))
c           ezz = ezz1*ezz2*zm2i(2)/den_z
           ezz = 0d0

           maa = mvv/den_a
           mzz = mvv/den_z

           !photon exchange:
           mat(k,isig1,isig3,j,l,1) = 
     4    	  maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)

           !Z exchange:
           mat(k,isig1,isig3,j,l,2) = 
     4          (mzz-ezz)*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)

	   endif !jlog
	 
           if (jlog2) then
              j = 2  ! never have g from upper line	   
              mat(k,isig1,isig3,j,l,1:2) = 0d0
	   endif !jlog
c
          enddo !isig
        enddo !isig
      enddo !l

c -----------------------
        
      else !cc_type
c
c charged current (k=5,6):
         den_w = dcmplx(q12(4,3)-xm2(3),xmg(3))
c      
      do l = 1,2   ! gluon polarization (kart. basis)
      
      if (jlog3) then
      j =  3 ! g from lower line
            
      if (k.eq.3) then 
         mww = dotcc(jqq(0,-1,1,is1,0),jqq(0,-1,2,is3,l))/den_w
         mat(k,-1,-1,j,l,3) = mww*clr(ifl(1,k),3,-1)*clr(ifl(3,k),3,-1)
      endif !k

      endif !jlog2
  
      if (jlog2) then
         j =  2                 ! no g from upper line 
         mat(k,-1,-1,j,l,3) = 0d0
      endif !jlog3
      enddo ! l-loop

      endif !nc/cc

c-----------------------------------------------------------------------      
c-----------------------------------------------------------------------
c ...IB: write the gluon helicity eigenstates in terms of the cartesian
c ...    polarizations ---> needed to calculate polarized xsection

      if (ipol.eq.1) then
       mat_ghe(:,:,:,:,0,:) = mat(:,:,:,:,0,:)

       mat_ghe(:,:,:,:,1,:) = dcmplx(-1./dsqrt(2.d0),0.d0)
     &                      *(mat(:,:,:,:,1,:)
     &                      -dcmplx(0.d0,gsign*1.d0)
     &                      *mat(:,:,:,:,2,:))

       mat_ghe(:,:,:,:,2,:) = dcmplx(1./dsqrt(2.d0),0.d0)
     &                      *(mat(:,:,:,:,1,:)
     &                      +dcmplx(0.d0,gsign*1.d0)
     &                      *mat(:,:,:,:,2,:)) 
c
       mat(:,:,:,:,:,:) = mat_ghe(:,:,:,:,:,:)
      endif
c ----------------------------------------------------------------------      
c-----------------------------------------------------------------------
c          
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
     
 999   continue

        if(zonly) mat(:,:,:,:,:,1) = 0d0 ! no gamma exchange
        if(aonly) mat(:,:,:,:,:,2) = 0d0 ! no Z exchange
c
c        if(zonly) mat_ghe(:,:,:,:,:,1) = 0d0 ! no gamma exchange
c        if(aonly) mat_ghe(:,:,:,:,:,2) = 0d0 ! no Z exchange        

c sum the graphs, square them and map them onto uucc, uuss etc.
      if (nc_type) then
         res = 0d0
         mm = 0
 	  do j = 2,3     
        do isig1 = -1,1,2
           do isig3 = -1,1,2
 	    	  do l = 1,2
              	     mm(k,isig1,isig3,j,l) = 0
               	     do i = 1,3
                        mm(k,isig1,isig3,j,l) = 
     1                     mm(k,isig1,isig3,j,l) + 
     1                     (mat(k,isig1,isig3,j,l,i))
c
c ...IB: switch to helicity eigenstates basis
c                        mm(k,isig1,isig3,j,l) = 
c     1                       mm(k,isig1,isig3,j,l) + 
c     1                     (mat_ghe(k,isig1,isig3,j,l,i))

                 enddo !i

c ...IB modified to account for polarized processes
                  if(ipol.eq.0) then                 
                     res(k,j) = res(k,j) 
     &                     + dreal(mm(k,isig1,isig3,j,l))**2
     &                     + dimag(mm(k,isig1,isig3,j,l))**2
                    elseif(ipol.eq.1) then
                       if(gsign.eq.1) then 
c
c ...IB  Initial quark
c ...    (note that for antiparticles isig = -physical helicity
c ...    --> added extra factor is1*is3)
c                    
                       res(k,j) = res(k,j) + sign(1)*sign(3)*( 
     &                      isig1*isig3*dreal(mm(k,isig1,isig3,j,l))**2
     &                    + isig1*isig3*dimag(mm(k,isig1,isig3,j,l))**2)
c                       
                       else
c ...IB  Initial gluon
c ...    (need to re-write the cartesian polarization in terms of 
c...      helicity eigenstates)
c ...    (added extra factor is1 for initial antiparticle)
c
                      res(k,j) = res(k,j) + sign(1)*( 
     &                  isig1*(-2*l+3)*dreal(mm(k,isig1,(isig3),j,l))**2
     &                 +isig1*(-2*l+3)*dimag(mm(k,isig1,isig3,j,l))**2)
c                     
                       endif !inital parton
c                           
                     else 
                        print *, "Polarization can be only 0 or 1,",
     $                       " aborting"
                        STOP
                     endif   !polarization
              enddo !l
          enddo !isig3		     
        enddo !isig1
c         res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
         res(k,j) = res(k,j)*4d0*fpials(j)   ! C_2*9 is the color factor
      enddo !j
c --------------------------------------------------------------------------
      else ! cc
      
        res = 0d0
        mm = 0d0
 	   do j = 2,3     
 	    do l = 1,2
               do i = 1,3
                  mm(k,-1,-1,j,l) = 
     1                mm(k,-1,-1,j,l) + 
     1                (mat(k,-1,-1,j,l,i))
c                  mm(k,-1,-1,j,l) = 
c     1                mm(k,-1,-1,j,l) + 
c     1                (mat_ghe(k,-1,-1,j,l,i))
               enddo !i
c
c ...IB modified to account for polarized processes
                  if(ipol.eq.0) then 

                  res(k,j) = res(k,j) 
     &               + dreal(mm(k,-1,-1,j,l))**2
     &               + dimag(mm(k,-1,-1,j,l))**2

                  elseif(ipol.eq.1) then
                    if(gsign.eq.1) then
c
c ...IB  Initial quark
                    res(k,j) = res(k,j) + sign(1)*sign(3)*(
     &                 dreal(mm(k,-1,-1,j,l))**2
     &               + dimag(mm(k,-1,-1,j,l))**2)
c
c ...IB  Initial gluon
                    else
                    res(k,j) = res(k,j) + sign(1)*(-1*(-2*l+3))*(
     &                 dreal(mm(k,-1,-1,j,l))**2
     &               + dimag(mm(k,-1,-1,j,l))**2)
c
                    endif !inital parton
c                           
                   else 
                    print*, 'set polarization in powheg.input'
                   endif   !polarization
c
        enddo !l
c        res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
        res(k,j) = res(k,j)*4d0*fpials(j)   ! C_2*9 is the color factor
       enddo !j
                  
      endif !nc/cc

c
c  set all processes to zero if photon virtuality falls below cutoff
      if (qdamp) then 
      if (abs(q12(4,3)).lt.qsqAmin) then       
         res(k,3) = res(k,3)*0.5d-20
         res(k,2) = res(k,2)*0.5d-20
      endif
      endif

            
      if (jmin.eq.3) then
         ans(2) = 0d0
      elseif (jmax.eq.2) then
         ans(3) = 0d0
      endif

      do j = 2,3
         ans(j) = res(k,j)
      enddo
      ans(1) = res(k,3)+res(k,2)

      return
      end

