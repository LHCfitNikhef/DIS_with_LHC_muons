      subroutine mboost(m,vec,beta,vin,vout)
c     boosts the m vectors vin(0:3,m) into the vectors vout(0:3,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(0:3,m),vout(0:3,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      real * 8 tiny
      parameter (tiny=1d-14)
c      if (abs(beta).ge.1d0) then
c         write(*,*) '********** WARNING ***************'
c         write(*,*) 'mboost called with beta=',beta
c         write(*,*) '**********************************'
c     endif
      if(beta == 0) then
         vout=vin
         return
      endif
      if (beta.ge.1d0) then
         beta = 1-tiny
      elseif (beta.le.-1d0) then
         beta = -1+tiny
      endif
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart) = vin(idim,ipart)-vec(idim)*vdotb
     $           +vec(idim)*gamma*(vdotb+beta*vin(0,ipart))
         enddo
         vout(0,ipart)=gamma*(vin(0,ipart)+vdotb*beta)
      enddo
      end

      subroutine mrotate(dir,sinphi,cosphi,vec)
c Rotates vector vec counterclockwise around the direction
c dir (|dir|=1) with angle phi, given sin phi and cos phi.
      implicit none
      real * 8 sinphi,cosphi,dir(3),vec(3)
      real * 8 dircrossvec(3),dirdotvec
      integer i
      dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
      dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
      dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
      dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
      do i=1,3
         vec(i)=vec(i)+sinphi*dircrossvec(i)
     #        -(1-cosphi)*(vec(i)-dir(i)*dirdotvec)
      enddo
      end

      subroutine rotate3tovec(v,vec)
c Rotate vec with the rotation that brings the positive
c third direction along v.
      implicit none
      real * 8 v(3),vec(3)
      real * 8 d(3),dir(3),cosphi,sinphi
      d = v/sqrt(v(1)**2+v(2)**2+v(3)**2)
      cosphi = d(3)
      dir(1) = d(2)
      dir(2) = -d(1)
      dir(3) = 0
      sinphi = sqrt(d(1)**2+d(2)**2)
      if(sinphi.lt.1d-9) return
      dir = - dir/sinphi
      call mrotate(dir,sinphi,cosphi,vec)
      end

      subroutine boost2reson(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(0:3),pin(0:3,nm),pout(0:3,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(0)
      if(beta < 1d-30) then
         pout = pin
      else
         vec(1)=pres(1)/(beta*pres(0))
         vec(2)=pres(2)/(beta*pres(0))
         vec(3)=pres(3)/(beta*pres(0))
         call mboost(nm,vec,-beta,pin,pout)
      endif
      end

      subroutine boost2resoninv(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(0:3),pin(0:3,nm),pout(0:3,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(0)
      vec(1)=pres(1)/(beta*pres(0))
      vec(2)=pres(2)/(beta*pres(0))
      vec(3)=pres(3)/(beta*pres(0))
      call mboost(nm,vec,beta,pin,pout)
      end

      !     Q in the lab frame (normal DIS variable), p in the lab frame and
!     pout in the breit frame. This uses appendix seven of the DIS book
!     by A Cooper-Sakar. Page 208 Appendix 7.11. 
      subroutine mlab2breit2(m,Q,p,pout,isPlus)
      implicit none
      integer i,j,m
      logical isPlus
      double precision Q(0:3), p(0:3,m), pin(0:3,m), pout(0:3,m), Qval
      double precision cosphi, sinphi, invQval, norm, q0, q1, q2, q3
      double precision z(3), Qb(0:3), modrot
      parameter (z = (/0d0,0d0,1d0/))
      double precision invmsq
      double precision trans(0:3,0:3)

      Qval = DSqrt(-invmsq(Q))
      invQval = 1d0/Qval
      Qb = Q                    ! Local copy
      pin = p                   ! Local copy
!     The routine assumes "+" like incoming direction. If that is not
!     the case, we revert the z-component.
      q0 = Qb(0)
      q1 = Qb(1)
      q2 = Qb(2)
      q3 = Qb(3)
      if(.not.isPlus) then
         q3 = -Qb(3)
      endif
      norm = 1d0/(q0 - q3)


!     Compute composite Lorentz transformation
      trans(0,0) = q0*invQval + Qval*norm
      trans(0,1) = -q1*invQval
      trans(0,2) = -q2*invQval
      trans(0,3) = -q3*invQval -Qval*norm

      trans(1,0) = -q1*norm
      trans(1,1) = 1d0
      trans(1,2) = 0d0
      trans(1,3) = q1*norm

      trans(2,0) = -q2*norm
      trans(2,1) = 0d0
      trans(2,2) = 1d0
      trans(2,3) = q2*norm

      trans(3,0) = q0*invQval
      trans(3,1) = -q1*invQval
      trans(3,2) = -q2*invQval
      trans(3,3) = -q3*invQval
      
      if (.not.isPlus) then
         trans(:,3) = - trans(:,3)
         trans(3,:) = - trans(3,:)
      endif
      
!     Perform transformation
      do j=1,m
         pout(:,j) = MATMUL(trans,pin(:,j))
      enddo
      end

!     Q in the lab frame (normal DIS variable), p in the breit frame and
!     pout in the lab frame. This uses appendix seven of the DIS book by
!     A Cooper-Sakar. Page 208 Appendix 7.11. This is the inverse of the
!     transformation mlab2breit2.
      subroutine mbreit2lab2(m,Q,p,pout,isPlus)
      implicit none
      integer i,j,m
      logical isPlus
      double precision Q(0:3), p(0:3,m), pin(0:3,m), pout(0:3,m), Qval
      double precision cosphi, sinphi, norm, invQval
      double precision z(3), Qb(0:3), modrot, q0, q1, q2, q3
      parameter (z = (/0d0,0d0,1d0/))
      double precision invmsq
      double precision trans(0:3,0:3)

      Qval = DSqrt(-invmsq(Q))
      invQval = 1d0/Qval
      Qb = Q                    ! Local copy
      pin = p                   ! Local copy
!     The routine assumes "+" like incoming direction. If that is not
!     the case, we revert the z-component.
      q0 = Qb(0)
      q1 = Qb(1)
      q2 = Qb(2)
      q3 = Qb(3)
      if(.not.isPlus) then
         q3 = -Qb(3)
      endif
      norm = 1d0/(q0 - q3)


!     Compute composite Lorentz transformation
      trans(0,0) = q0*invQval + Qval*norm
      trans(0,1) = q1*norm
      trans(0,2) = q2*norm
      trans(0,3) = -q0*invQval

      trans(1,0) = q1*invQval
      trans(1,1) = 1d0
      trans(1,2) = 0d0
      trans(1,3) = -q1*invQval

      trans(2,0) = q2*invQval
      trans(2,1) = 0d0
      trans(2,2) = 1d0
      trans(2,3) = -q2*invQval

      trans(3,0) = Qval*norm + q3*invQval
      trans(3,1) = q1*norm
      trans(3,2) = q2*norm
      trans(3,3) = -q3*invQval
      
      if (.not.isPlus) then
         trans(:,3) = - trans(:,3)
         trans(3,:) = - trans(3,:)
      endif
      
!     Perform transformation
      do j=1,m
         pout(:,j) = MATMUL(trans,pin(:,j))
      enddo
      end

      double precision function invmsq(p)
      double precision p(0:3)

      invmsq = p(0)**2 - p(1)**2 - p(2)**2 - p(3)**2
      end

