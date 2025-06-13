!=======================================
!     AUXILIARY EVENT SHAPE ROUTINES (FROM LIBRARY EVSHPLIB) 
      
      subroutine evs_tauzE_DIS(parr,tau_zE)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp 
      real(dp), intent(in)  :: parr(:,:)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp)          :: Qvec(4)
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      
      real(dp), intent(out) :: tau_zE
      integer           :: n, i, j
      real(dp) :: temp
      Qvec=zero
      Qvec(3)=one
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
         if (parr(3,i)*Qvec(3) > zero) then
            n = n+1
            pcur(:,n) = parr(:,i)
            Ecur = Ecur + parr(4,i)
         end if
      end do
      
!--- assume that photon is along 3-direction
      EoverhalfQ = Ecur / (half*abs(Qvec(3)) )
      
!------ tau_z ----
      if (n == 0) then
         thr = zero
      else
         thr = sum(pcur(3,:n)) / (half * Qvec(3))
      end if
      
      if (Ecur /= zero) tau_zE = one - thr / EoverhalfQ

      end subroutine evs_tauzE_DIS


      subroutine evs_tauzq_DIS(parr,Qvec,tau_zQ)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp 
      real(dp), intent(in)  :: parr(:,:)
      real(dp), intent(in) :: Qvec(4)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur, Q
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: tau_zQ
      integer           :: n, i, j
      real(dp) :: temp
      Q = sqrt(dot_product(Qvec(:3),Qvec(:3))-Qvec(4)*Qvec(4))
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
          if (parr(3,i)*Qvec(3) > zero) then
             n = n+1
             pcur(:,n) = parr(:,i)
          end if
       end do
!------ tau_z ----
       if (n == 0) then
          thr = zero
       else
          thr = sum(pcur(3,:n)) / (half * Q)
       end if
       
       tau_zQ = one - thr

       end subroutine evs_tauzq_DIS      

 !======================================================================
      subroutine es_thr_tQ_selctd(pcur, Ecur, Qvec, n, thr, a)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp
      real(dp), intent(in)  :: pcur(:,:), Ecur, Qvec(:)
      integer,  intent(in)  :: n
      real(dp), intent(out) :: thr, a(3)
!------------------------------------------------
      real(dp) :: thrtmp, atmp(3),mod3sq
      integer  :: i, j, maxiter
      logical :: first_time = .true. 
      
      thr = zero
      select case (n)
      case (0)
!-- give a value for safety
         a=zero
         a(3) = Qvec(3)/abs(Qvec(3)) 
      case (1)
         thr = pcur(4,1)
         a   = pcur(:3,1) / thr
      case default
         maxiter = 2**(n-1) - 1
         if (first_time .and. n> 3) then 
            write(0,*) 'VERY INEFFICIENT ALGORTHM FOR tau_zQ 
     Cfor large n       !' 
            first_time  = .false.
         end if
!-- algorithm tries all sign combinations of the momenta (except for
! overall sign which is adjusted at the end) so as to get all potential
! thrust axes, and then calculates the thrust for each one.
!     
         thr = zero
         do i = 0, maxiter
            atmp = pcur(:3,1)
            do j = 2, n
               if (.not. btest(i,j-2)) then
                  atmp = atmp + pcur(:3,j)
               else
                  atmp = atmp - pcur(:3,j)
               end if
            end do
!-- on 10^4 events, this seems to give identical results to below.
            thrtmp = mod3sq(atmp)
            if (thrtmp > thr) then
               thr = thrtmp
               a   = atmp
            end if
         end do
!-- take square roots, normalise, and get axis in right direction
!-- thr should always be > 0?
         if (thr <= zero) then
            write(0,*) 'WARNING: thr^2 had illegal value:',(thr)
            thr = zero
            a = (/zero,zero,one/)
         else
            thr = sqrt(thr)
            a = a / thr
         end if
         
       if (dot_product(a,Qvec(:3)) < 0)  a = a*(-one)
       
      end select
      thr = thr / (half*abs(Qvec(3)))
      end subroutine es_thr_tQ_selctd
      
      

      subroutine evs_tautE_DIS(parr,tau_tE)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp
      real(dp), intent(in)  :: parr(:,:)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp)          :: Qvec(4)
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: tau_tE
      integer           :: n, i, j
      real(dp) :: temp

      
      interface
      subroutine es_thr_tQ_selctd(pcur, Ecur, Qvec, n, thr, a)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp), intent(in)  :: pcur(:,:), Ecur, Qvec(:)
      integer,  intent(in)  :: n
      real(dp), intent(out) :: thr, a(3)
      end       
      end interface 

      
      Qvec=zero
      Qvec(3)=one
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
         if (parr(3,i)*Qvec(3) > zero) then
            n = n+1
            pcur(:,n) = parr(:,i)
            Ecur = Ecur + parr(4,i)
         end if
      end do

!--- assume that photon is along 3-direction
      EoverhalfQ = Ecur / (half*abs(Qvec(3)) )
      
!------ tau_t ----
      call es_thr_tQ_selctd(pcur, Ecur, Qvec, n, thr, a)
      if (Ecur /= zero) tau_tE = one - thr / EoverhalfQ
      
      end subroutine evs_tautE_DIS      

! squared modulus of 3 vector
      function mod3sq(p) result(sqr)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), intent(in) :: p(1:3)
      real(dp) :: sqr
      sqr = p(1)**2 + p(2)**2 + p(3)**2
      end function mod3sq
      

      subroutine evs_brdzE_DIS(parr, brd_zE)
      implicit none
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp
      real(dp), intent(in)  :: parr(:,:)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp)          :: Qvec(4)
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: brd_zE
      integer           :: n, i, j
      real(dp) :: temp

      Qvec=zero
      Qvec(3)=one
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
          if (parr(3,i)*Qvec(3) > zero) then
             n = n+1
             pcur(:,n) = parr(:,i)
             Ecur = Ecur + parr(4,i)
          end if
       end do
       
!------ brd_zE ---
       select case(n)
      case(0)
         brd = zero             !es_undefined_value
      case default
         brd = zero
         do i = 1, n
!brd = brd + sqrt(mod3sq(pcur(:3,i) - a*dot_product(a,pcur(:3,i))))
!--   assume z axis is along 3 direction!
            brd = brd + sqrt(pcur(1,i)**2 + pcur(2,i)**2)
         end do
         brd = half*brd / Ecur
         brd_zE = brd
      end select    
      
      
      end subroutine evs_brdzE_DIS
      
      
      subroutine evs_rhoE_DIS(parr, rho_E)
      implicit none
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp), parameter :: two = 2.0_dp, three = 3.0_dp
      real(dp), intent(in)  :: parr(:,:)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp)          :: Qvec(4)
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: rho_E
      integer           :: n, i, j
      real(dp) :: temp,dot_2v
      
      Qvec=zero
      Qvec(3)=one
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
         if (parr(3,i)*Qvec(3) > zero) then
            n = n+1
            pcur(:,n) = parr(:,i)
            Ecur = Ecur + parr(4,i)
         end if
      end do
      
!------rho_E ----
      select case(n)
      case(0)
         rho = zero             !es_undefined_value
      case default
!     ptot = sum(pcur(:4,:n),dim=2)
         do i=1,4
            ptot(i) = sum(pcur(i,:n))
         end do
         rho  = max(zero,dot_2v(ptot,ptot))
!     write(0,'(4f15.8)') real(pcur(:4,:n))
         rho = rho / (two*Ecur)**2
      end select
      rho_E = rho
      
      
      end subroutine evs_rhoE_DIS
      

      subroutine evs_CprE_DIS(parr,Cpr_E)
      implicit none
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp), parameter :: two = 2.0_dp, three = 3.0_dp
      real(dp), intent(in)  :: parr(:,:)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp)          :: Qvec(4)
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: Cpr_E
      integer           :: n, i, j
      real(dp) :: temp
      
      Qvec=zero
      Qvec(3)=one
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
         if (parr(3,i)*Qvec(3) > zero) then
            n = n+1
            pcur(:,n) = parr(:,i)
            Ecur = Ecur + parr(4,i)
         end if
      end do
      
!------ cpr_E ----
      select case(n)
      case(0)
         Cpr = zero             !es_undefined_value
      case(1)
         Cpr = zero
      case default
         Cpr = zero
         do i = 1, n-1
            do j = i+1, n
               costheta = dot_product(pcur(:3,i),pcur(:3,j))/
     C              (pcur(4,i)*pcur(4,j))
               Cpr = Cpr + (pcur(4,i)*pcur(4,j))*(one - costheta**2)
            end do
         end do
!--- factor of half cancelled by not doing symmetrized sum
         Cpr = three * Cpr / (Ecur**2)
      end select
      cpr_E = cpr
      
      end subroutine evs_CprE_DIS
      
      function dot_2v(p,q) result(dot)
      implicit none
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), intent(in) :: p(1:4),q(1:4)
      real(dp) :: dot
      dot = p(4)*q(4) - p(3)*q(3) - p(2)*q(2) - p(1)*q(1)
      end function dot_2v
