MODULE hist2d

   !TODO 
   ! add underflow and overflow bins to the histogram structure.

   IMPLICIT NONE

   !PRIVATE
   !PUBLIC :: filldd,new_ddDist,lineq_ddDist,logeq_ddDist,bookupeqlogdd

   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   ! Contains the histogram for the double differential distribution.
   TYPE ddDist
      CHARACTER(LEN=80) :: str
      INTEGER(KIND=4) :: nbinx,nbiny
      REAL(KIND=dp),DIMENSION(:,:,:),ALLOCATABLE :: hist
      REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: xedges,yedges
      REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: xwidth,ywidth
      INTEGER(KIND=4) :: nhits
      REAL(KIND=dp) :: nevents

      CONTAINS
         PROCEDURE :: fill
         PROCEDURE :: is_present
         PROCEDURE :: write_to
   END TYPE ddDist

   ! All the booked histograms are stored in dists.
   TYPE(ddDist),DIMENSION(:),POINTER,SAVE :: dists

   INTERFACE check_histogram_edges
      MODULE PROCEDURE check_histogram_edges_realdp
      MODULE PROCEDURE check_histogram_edges_int4
   END INTERFACE check_histogram_edges

   CONTAINS

      SUBROUTINE init_dists
         INTEGER(KIND=4) :: i
         ALLOCATE(dists(1))
         DO i=LBOUND(dists,DIM=1),UBOUND(dists,DIM=1)
            dists(i)%str=""
         END DO
      END SUBROUTINE init_dists

      SUBROUTINE ddDist_add_event
         INTEGER(KIND=4) :: i
         DO i=LBOUND(dists,DIM=1),UBOUND(dists,DIM=1)
            dists(i)%nevents=dists(i)%nevents+1.0_dp
         END DO
      END SUBROUTINE ddDist_add_event


      SUBROUTINE filldd(str,sig,x,y)
         CHARACTER(LEN=*),INTENT(IN) :: str
         REAL(KIND=dp),INTENT(IN) :: x,y
         REAL(KIND=dp),DIMENSION(*),INTENT(IN) :: sig
         INTEGER(KIND=4) :: i

         WRITE(*,*) "DEBUG filldd: ",str
         i=get_dist_index(str)
         CALL dists(i)%fill(sig,x,y)
      END SUBROUTINE filldd
      
      INTEGER(KIND=4) FUNCTION get_dist_index(str) RESULT(ind)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4) :: i
         WRITE(*,*) "DEBUG get_dist_index:",str
         DO i=LBOUND(dists,DIM=1),UBOUND(dists,DIM=1)
            IF (str.EQ.TRIM(ADJUSTL(dists(i)%str)))THEN
               ind=i
               RETURN
            END IF
         END DO
      END FUNCTION get_dist_index

      SUBROUTINE bookupeqlindd(str,lowx,upx,nbinx,lowy,upy,&
                                       & nbiny,nwgt)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: nbinx,nbiny,nwgt
         REAL(KIND=dp),INTENT(IN) :: lowx,upx,lowy,upy
         INTEGER(KIND=4) :: i
         LOGICAL :: booked
         TYPE(ddDist),DIMENSION(:),POINTER :: tmp

         booked=.FALSE.
         DO i=1,UBOUND(dists,DIM=1)
            IF (dists(i)%is_present(str).EQ.1) THEN
               dists(i)=lineq_ddDist(str,lowx,upx,nbinx,lowy,upy,nbiny,nwgt)
               booked=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.booked) THEN
            WRITE(*,*) "Allocating new distribution"
            ALLOCATE(tmp(UBOUND(dists,DIM=1)))
            tmp=dists
            DEALLOCATE(dists)
            ALLOCATE(dists(UBOUND(tmp,DIM=1)+1))
            dists(1:UBOUND(tmp,DIM=1))=tmp(:)
            DEALLOCATE(tmp)
            dists(UBOUND(dists,DIM=1))=lineq_ddDist(str,lowx,upx,nbinx,lowy,upy,nbiny,nwgt)
         END IF

      END SUBROUTINE bookupeqlindd

      SUBROUTINE bookupeqlogdd(str,lowx,upx,nbinx,lowy,upy,&
                                       & nbiny,nwgt)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: nbinx,nbiny,lowx,upx,lowy,upy,nwgt
         INTEGER(KIND=4) :: i
         LOGICAL :: booked
         TYPE(ddDist),DIMENSION(:),POINTER :: tmp

         booked=.FALSE.
         DO i=1,UBOUND(dists,DIM=1)
            IF (dists(i)%is_present(str).EQ.1) THEN
               dists(i)=logeq_ddDist(str,lowx,upx,nbinx,lowy,upy,nbiny,nwgt)
               booked=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.booked) THEN
            WRITE(*,*) "Allocating new distribution"
            ALLOCATE(tmp(UBOUND(dists,DIM=1)))
            tmp=dists
            DEALLOCATE(dists)
            ALLOCATE(dists(UBOUND(tmp,DIM=1)+1))
            dists(1:UBOUND(tmp,DIM=1))=tmp(:)
            DEALLOCATE(tmp)
            dists(UBOUND(dists,DIM=1))=logeq_ddDist(str,lowx,upx,nbinx,lowy,upy,nbiny,nwgt)
         END IF

      END SUBROUTINE bookupeqlogdd

      SUBROUTINE bookupeqbpdlogdd(str,lowx,upx,nbpdx,lowy,upy,&
                                       & nbpdy,nwgt)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: nbpdx,nbpdy,lowx,upx,lowy,upy,nwgt
         INTEGER(KIND=4) :: i
         LOGICAL :: booked
         TYPE(ddDist),DIMENSION(:),POINTER :: tmp

         booked=.FALSE.
         DO i=1,UBOUND(dists,DIM=1)
            IF (dists(i)%is_present(str).EQ.1) THEN
               dists(i)=logeqperdecade_ddDist(str,lowx,upx,nbpdx,lowy,&
                                                & upy,nbpdy,nwgt)
               booked=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.booked) THEN
            WRITE(*,*) "Allocating new distribution"
            ALLOCATE(tmp(UBOUND(dists,DIM=1)))
            tmp=dists
            DEALLOCATE(dists)
            ALLOCATE(dists(UBOUND(tmp,DIM=1)+1))
            dists(1:UBOUND(tmp,DIM=1))=tmp(:)
            DEALLOCATE(tmp)
            dists(UBOUND(dists,DIM=1))=logeqperdecade_ddDist(str,lowx,&
                                                & upx,nbpdx,lowy,&
                                                & upy,nbpdy,nwgt)
         END IF

      END SUBROUTINE bookupeqbpdlogdd

      INTEGER FUNCTION is_present(self,str)
         CLASS(ddDist),INTENT(IN) :: self
         CHARACTER(LEN=80),INTENT(IN) :: str
         IF(TRIM(self%str).EQ."".OR.self%str.EQ.str) THEN
            is_present=1
         ELSE
            is_present=2
         END IF
      END FUNCTION is_present

      ! Book a generic 2D histogram providing arrays for the edges.
      TYPE(ddDist) FUNCTION new_ddDist(str,xedges,yedges,nwgt) RESULT(dist)
         CHARACTER(LEN=*),INTENT(IN) :: str
         REAL(KIND=dp),DIMENSION(:),INTENT(IN) :: xedges,yedges
         INTEGER(KIND=4),INTENT(IN) ::nwgt
         REAL(KIND=dp) :: inc
         INTEGER(KIND=4) :: i

         dist%str=TRIM(ADJUSTL(str))
         dist%nhits=0
         dist%nevents=0.0_dp
         CALL check_histogram_edges(str,xedges)
         CALL check_histogram_edges(str,yedges)
         ALLOCATE(dist%xedges,SOURCE=xedges)
         ALLOCATE(dist%yedges,SOURCE=yedges)
         ALLOCATE(dist%xwidth(SIZE(xedges)-1))
         ALLOCATE(dist%ywidth(SIZE(yedges)-1))
         CALL binwidth_from_edges(str,dist%xedges,dist%xwidth)
         CALL binwidth_from_edges(str,dist%yedges,dist%ywidth)
         dist%nbinx=SIZE(dist%xwidth)
         dist%nbiny=SIZE(dist%ywidth)
         ALLOCATE(dist%hist(0:dist%nbinx+1,0:dist%nbiny+1,nwgt))
         dist%hist=0.0_dp
      END FUNCTION new_ddDist

      ! Book a 2D histogram with equal spaced bins in linear space.
      TYPE(ddDist) FUNCTION lineq_ddDist(str,lowx,upx,nbinx,lowy,upy,&
                                       & nbiny,nwgt) RESULT(dist)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: nbinx,nbiny,nwgt
         REAL(KIND=dp),INTENT(IN) :: lowx,upx,lowy,upy
         REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: xedges,yedges
         
         CALL check_histogram_edges(str,(/lowx,upx/))
         CALL check_histogram_edges(str,(/lowy,upy/))
         ALLOCATE(xedges(nbinx+1))
         ALLOCATE(yedges(nbinx+1))
         CALL compute_edges(lowx,upx,xedges)
         CALL compute_edges(lowy,upy,yedges)
         dist=new_ddDist(str,xedges,yedges,nwgt)
         DEALLOCATE(xedges,yedges)

         CONTAINS
            SUBROUTINE compute_edges(low,up,edges)
               REAL(KIND=dp),INTENT(IN) :: low,up
               REAL(KIND=dp),DIMENSION(:),INTENT(OUT) :: edges
               REAL(KIND=dp) :: inc
               INTEGER(KIND=4) :: i
               inc=(REAL(up,KIND=dp)-REAL(low,KIND=dp)) &
                 & /(REAL(SIZE(edges,1),KIND=dp)-1.0_dp)
               DO i=LBOUND(edges,1),UBOUND(edges,1)
                  edges(i)=low+(i-1)*inc
               END DO
            END SUBROUTINE compute_edges

      END FUNCTION lineq_ddDist

      ! Book a 2D histogram with equal spaced bins in logspace
      TYPE(ddDist) FUNCTION logeq_ddDist(str,lowx,upx,nbinx,lowy,upy,nbiny,nwgt) RESULT(dist)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: lowx,upx,nbinx,lowy,upy,nbiny,nwgt
         REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: xedges,yedges

         CALL check_histogram_edges(str,(/lowx,upx/))
         CALL check_histogram_edges(str,(/lowy,upy/))
         ALLOCATE(xedges(nbinx+1))
         ALLOCATE(yedges(nbiny+1))
         CALL compute_edges_log(lowx,upx,xedges)
         CALL compute_edges_log(lowy,upy,yedges)
         dist=new_ddDist(str,xedges,yedges,nwgt)
         DEALLOCATE(xedges,yedges)

         CONTAINS
            SUBROUTINE compute_edges_log(low,up,edges)
               IMPLICIT NONE
               INTEGER(KIND=4),INTENT(IN) :: low,up
               REAL(KIND=dp),DIMENSION(:),INTENT(INOUT) :: edges
               INTEGER(KIND=4) :: i
               REAL(KIND=dp) :: inc
               inc=(REAL(up,KIND=dp)-REAL(low,KIND=dp)) &
                 & /(REAL(SIZE(edges,1),KIND=dp)-1.0_dp)
               DO i=LBOUND(edges,1),UBOUND(edges,1)
                  edges(i)=10.0_dp**(low+(i-1)*inc)
               END DO
            END SUBROUTINE compute_edges_log
      END FUNCTION logeq_ddDist

      TYPE(ddDist) FUNCTION logeqperdecade_ddDist(str,lowx,upx,nbpdx,lowy,&
                                                & upy,nbpdy,nwgt) RESULT(dist)
         CHARACTER(LEN=*),INTENT(IN) :: str
         INTEGER(KIND=4),INTENT(IN) :: lowx,upx,nbpdx,lowy,upy,nbpdy,nwgt
         INTEGER(KIND=4) :: nbinx,nbiny
         REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: xedges,yedges

         CALL check_histogram_edges(str,(/lowx,upx/))
         CALL check_histogram_edges(str,(/lowy,upy/))
         !Add bins from bpd
         nbinx=(upx-lowx)*nbpdx
         nbiny=(upy-lowy)*nbpdy
         dist%nbinx=nbinx
         dist%nbiny=nbiny
         ALLOCATE(xedges(dist%nbinx+1))
         ALLOCATE(yedges(dist%nbiny+1))
         CALL compute_edges_logpd(lowx,upx,nbpdx,xedges)
         CALL compute_edges_logpd(lowy,upy,nbpdy,yedges)
         dist=new_ddDist(str,xedges,yedges,nwgt)
         DEALLOCATE(xedges,yedges)

         CONTAINS
            SUBROUTINE compute_edges_logpd(low,up,bpd,edges)
               INTEGER(KIND=4),INTENT(IN) :: low,up,bpd
               REAL(KIND=dp),DIMENSION(:),INTENT(OUT) :: edges
               INTEGER(KIND=4) :: i
               REAL(KIND=dp) :: inc
               inc=1.0_dp/REAL(bpd,KIND=dp)
               DO i=LBOUND(edges,1),UBOUND(edges,1)
                  edges(i)=10**(low+(i-1)*inc)
               END DO
            END SUBROUTINE compute_edges_logpd
      END FUNCTION logeqperdecade_ddDist

      SUBROUTINE fill(self,xsec,x,y)
         CLASS(ddDist),INTENT(INOUT) :: self
         REAL(KIND=dp),DIMENSION(SIZE(self%hist,3)),INTENT(IN) :: xsec
         REAL(KIND=dp),INTENT(IN) :: x,y
         REAL(KIND=dp) :: dxdy
         INTEGER(KIND=4) i,j

         i=findbin(x,self%xedges)
         j=findbin(y,self%yedges)

         WRITE(*,*) "DEBUG: Writing value to histogram "//TRIM(self%str)
         WRITE(*,*) "DEBUG: indices ",i,j
         WRITE(*,*) "DEBUG: upper bounds ",UBOUND(self%hist,1),UBOUND(self%hist,2)
         WRITE(*,*) "DEBUG: is allocated",ALLOCATED(self%hist)
         IF(i.EQ.0)THEN ! Underflow in x
            dxdy=self%xwidth(i)*self%ywidth(j)
            self%hist(i,j,:)=self%hist(i,j,:)+xsec/dxdy
            RETURN
         ELSE IF(j.EQ.0)THEN  ! Underflow in y
            dxdy=self%xwidth(i)*self%ywidth(j)
            self%hist(i,j,:)=self%hist(i,j,:)+xsec/dxdy
            RETURN
         ELSE IF(i.EQ.UBOUND(self%hist,1))THEN ! Overflow in x
            dxdy=self%xwidth(i)*self%ywidth(j)
            self%hist(i,j,:)=self%hist(i,j,:)+xsec/dxdy
            RETURN
         ELSE IF(j.EQ.UBOUND(self%hist,2))THEN ! Overflow in y
            dxdy=self%xwidth(i)*self%ywidth(j)
            self%hist(i,j,:)=self%hist(i,j,:)+xsec/dxdy
            RETURN
         ELSE ! Hit in x and y
            dxdy=self%xwidth(i)*self%ywidth(j)
            self%hist(i,j,:)=self%hist(i,j,:)+xsec/dxdy
            self%nhits=self%nhits+1
            RETURN
         END IF

         CONTAINS 

            INTEGER(KIND=4) FUNCTION findbin(val,arr) RESULT(bin)
               REAL(KIND=8),INTENT(IN) :: val
               REAL(KIND=8),DIMENSION(:),INTENT(IN) :: arr
               INTEGER(KIND=4) :: i
               bin=UBOUND(arr,DIM=1)
               DO i=LBOUND(arr,DIM=1),UBOUND(arr,DIM=1)
                  IF(val.LT.arr(i))THEN
                     bin=i-1
                     EXIT
                  END IF
               END DO
            END FUNCTION  findbin

      END SUBROUTINE fill


      SUBROUTINE check_histogram_edges_realdp(str,edges)
         CHARACTER(LEN=80),INTENT(IN) :: str
         REAL(KIND=dp),DIMENSION(:),INTENT(IN) :: edges
         CHARACTER(LEN=80) :: msg
         INTEGER(KIND=4) :: i
         REAL(KIND=dp) :: low,up
         DO i=LBOUND(edges,DIM=1),UBOUND(edges,DIM=1)-1
            low=edges(i)
            up=edges(i)
            IF(low.GT.up)THEN
               msg="ERROR in "//TRIM(ADJUSTL(str))// &
                  &": Lower limit larger than upper limit."
               ERROR STOP msg
            END IF
         END DO
      END SUBROUTINE check_histogram_edges_realdp

      SUBROUTINE check_histogram_edges_int4(str,edges)
         CHARACTER(LEN=80),INTENT(IN) :: str
         INTEGER,DIMENSION(:),INTENT(IN) :: edges
         CHARACTER(LEN=80) :: msg
         INTEGER(KIND=4) :: i,low,up
         DO i=LBOUND(edges,1),UBOUND(edges,1)-1
            low=edges(i)
            up=edges(i)
            IF(low.GT.up)THEN
               msg="ERROR in "//TRIM(ADJUSTL(str))// &
                  &": Lower limit larger than upper limit."
               ERROR STOP msg
            END IF
         END DO
      END SUBROUTINE check_histogram_edges_int4

      SUBROUTINE binwidth_from_edges(str,edges,binwidth)
         CHARACTER(LEN=80),INTENT(IN) :: str
         REAL(KIND=dp),DIMENSION(:),INTENT(IN) :: edges
         REAL(KIND=dp),DIMENSION(:),INTENT(OUT) :: binwidth
         CHARACTER(LEN=80) :: msg
         INTEGER(KIND=4) :: i

         !IF(SIZE(edges,1)-1.NE.SIZE(binwidth,1))THEN
         !   msg="ERROR in "//TRIM(ADJUSTL(str))// &
         !      &": array edges and binwidth not compatible."
         !   ERROR STOP msg
         !END IF

         DO i=LBOUND(binwidth,1),UBOUND(binwidth,1)
            binwidth(i)=edges(i+1)-edges(i)
         END DO
      END SUBROUTINE binwidth_from_edges

      SUBROUTINE write_to(self,filestr)
         CLASS(ddDist),INTENT(IN) :: self
         CHARACTER(LEN=80),INTENT(IN) :: filestr
         CHARACTER(LEN=80) :: filename
         INTEGER(KIND=4) :: un,i,j,k,wgt

         DO wgt=LBOUND(self%hist,3),UBOUND(self%hist,3)
            WRITE(filename,"(A,I0)") TRIM(ADJUSTL(filestr))//"-W",wgt
            OPEN(NEWUNIT=un,FILE=TRIM(filename),STATUS="UNKNOWN")
            WRITE(un,*) "Distribution",TRIM(self%str)
            WRITE(*,*) "Distribution",TRIM(self%str)
            WRITE(*,*) "DEBUG",SIZE(self%hist,1),SIZE(self%hist,2),SIZE(self%hist,3)

            ! Here we do not write the underflow and overflow bins.
            DO i=LBOUND(self%hist,1)+1,UBOUND(self%hist,1)-1
               DO j=LBOUND(self%hist,2)+1,UBOUND(self%hist,2)-1
                  WRITE(un,"(5E16.7)") self%xedges(i),self%xedges(i+1),&
                & self%yedges(j),self%yedges(j+1),&
                & self%hist(i,j,wgt)/self%nevents
               END DO
            END DO
            CLOSE(un)
         END DO

      END SUBROUTINE write_to

END MODULE hist2d
