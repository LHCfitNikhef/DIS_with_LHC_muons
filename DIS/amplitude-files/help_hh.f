c***************************************************************************
      function dotrr(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      double precision dotrr,p1(0:3),p2(0:3)

      dotrr = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end
      
c ------------------------------------------------------------

      double complex function dotcc(v1,v2)
      implicit none
      double complex v1(0:3), v2(0:3)
      dotcc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


      double complex function dotrc(v1,v2)
      implicit none
      double precision v1(0:3)
      double complex  v2(0:3)
      dotrc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


      double complex function dotqj(v2)
      implicit none
      double precision v1(0:3)
      double complex  v2(0:5)
      v1(0) = dreal(v2(4))
      v1(1) = dreal(v2(5))
      v1(2) = dimag(v2(5))
      v1(3) = dimag(v2(4))
      dotqj = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


      double complex function contract_Tjj(T,j1,j2)
      implicit none
c contract complex rank 2 tensor T^{mu,nu} with two complex vectors j1 and j2
      complex*16 T(0:3,0:3), j1(0:3), j2(0:3), resv(0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Tjj = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
      return
      end

      subroutine contract_T1j(T,jc,jout)
      implicit none
c contract first index of complex rank 2 tensor T^{mu,nu} with complex vector jc
c jout(mu) = T^{nu,mu} jc(nu) 
      complex*16 T(0:3,0:3), jc(0:3), jout(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(0,mu)*jc(0) - T(1,mu)*jc(1) 
     &            - T(2,mu)*jc(2) - T(3,mu)*jc(3)
      enddo
      end



      subroutine contract_T2j(T,jc,jout)
      implicit none
c contract second index of complex rank 2 tensor T^{mu,nu} with complex vector jc
c jout(mu) = T^{mu,nu} jc(nu) 
      complex*16 T(0:3,0:3), jc(0:3), jout(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(mu,0)*jc(0) - T(mu,1)*jc(1) 
     &            - T(mu,2)*jc(2) - T(mu,3)*jc(3)
      enddo
      end
      
C-------------------------------------------------------------------

      FUNCTION SC3(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), AUX(0:3,3)
      REAL*8  A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

c--------------------------------------------------------------       
c       
      subroutine revers(ans,rev)
C  
c	rev(ans) exchanges ans(2) and ans(3);
c	needed in call of mg routines with switched quark lines
c
c	in:   ans(3) ... 3 dim array
c	out:  rev(3) ... 3 dim array (=ans switched)	
c
C  
      IMPLICIT NONE
      REAL*8    ans(3),rev(3)

	
	rev(1) = ans(1)
	rev(2) = ans(3)
	rev(3) = ans(2)
	
      end	
