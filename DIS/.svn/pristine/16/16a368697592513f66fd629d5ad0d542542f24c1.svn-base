      subroutine sample_x_lepton(random, xmin, xmax, xlep, jac) 
      implicit none
! input
      real *8 random       ! Random number
      real *8 xmin, xmax   ! Lower and upper limit on xlepton (computed in Born_phsp)
!     output
      real *8 xlep         ! Incoming lepton x fraction
      real *8 jac          ! Jacobian from xlep to random
      
!     REPLACE THIS with the importance sampling of your choice
      xlep = xmin + random * (xmax-xmin)
!     jacobian = ∂xlep/∂random
      jac = (xmax-xmin)
      end


      subroutine pdf_lepton_beam(ih, mufact2, x, pdf )
      implicit none
      include 'pwhg_pdf.h'
!     input
      integer ih                !Index of the incoming beam
      real * 8 mufact2          !factorisation scale (squared)
      real * 8 x                !x fraction
!     output
      real *8 pdf(-pdf_nparton:pdf_nparton)   !Array of pdfs

      pdf = 1e-8                !Do not set light quark pdfs to 0 otherwise the checklims is broken. This is in any case not used

!     COMMENT -- ADD THE FUNCTIONAL FORM YOU PREFER FOR THE FLAVOURS YOU NEED          
      pdf(11)  = exp(-1d0/x**2)
      pdf(-11)  = exp(-1d0/x**2)
      pdf(12)  = 1d0
      pdf(-12) = 1d0
      end
