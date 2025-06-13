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
      include 'pwhg_pdf.h'
      real * 8 lam5
      integer iord,iset,maxsets
      common/cgenericpdf/lam5,iord,iset,maxsets

      integer(kind=4),intent(in) :: ih        !Index of the incoming beam
      real(kind=8),intent(in)    :: mufact2   !factorisation scale (squared)
      real(kind=8),intent(in)    :: x         !x fraction
      real(kind=8),dimension(-pdf_nparton:pdf_nparton),
     1 intent(out)               :: pdf

      integer(kind=4),save       :: LEPpdf=-1
      logical,save               :: ini=.TRUE.

      interface
         real(kind=8) function powheginput(stringa)
            character(len=*),dimension(*),intent(in) :: stringa
         end function powheginput
         logical(kind=c_bool) function generic_has_id(iset,id) 
     1   bind(C,name="generic_has_id_")
            use iso_c_binding, only: c_int, c_bool
            integer(kind=c_int),intent(in) :: iset,id
         end function generic_has_id
      end interface

      pdf = 1e-8                ! Do not set light quark pdfs to 0 otherwise
                                ! the checklims is broken. This is in any 
                                ! case not used
      if(ini)then
         leppdf=int(powheginput("#LEPpdf"),kind=4)
         if(leppdf.lt.0)then
            write(*,"(A)") 
     1           "ERROR: Keyword LEPpdf missing in input card."
            stop
         end if
         call genericpdfset(LEPpdf)
         if(.NOT.generic_has_id(iset,ih))then
            write(*,"(A23,I4,A15)") "ERROR: Neutrino flavour",ih, 
     1                            " not available."
            stop
         end if
         ini=.FALSE.
      end if

      call genericpdfset(LEPpdf)
      call xf_pdgid(iset,ih,x,mufact2,pdf(ih))
      pdf(ih)=pdf(ih)/x

      end
