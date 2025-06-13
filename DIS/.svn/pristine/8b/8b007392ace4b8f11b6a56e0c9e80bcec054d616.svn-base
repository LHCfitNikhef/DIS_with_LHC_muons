!     DIS phase space taken from nlojet++
!     dPhi = dx dOmega/32 pi^2
!          = dx dcth/16pi * dphi/2pi       [cth = 1-2y]
!          = dx 2dy/16 pi * dphi/2pi  
!     with Q^2 = x y S      
      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'

      real * 8 xborn(ndiminteg-3)
      real * 8 vec(3), beta
! DIS variables
      real * 8 x,Q2,s,y, xl
      real * 8 xljac, jac, Eh, El
!     Limits on Q2 and x
      real * 8 Q2min, Q2max, xmin, xmax, ymin, ymax,ymn,ymx, xminl, xmaxl,
     $     xlm, xlp
      integer mu,k,pow, ixb
      logical ini
      data ini/.true./
      save ini, s, El, Eh, Q2min, Q2max, xmin, xmax, ymin, ymax, xminl, xmaxl

      real*8 powheginput,sqrts,theta,phi,cost,sint
      real*8 yaxis(3),zaxis(3),q(0:3)
      parameter (yaxis = (/0d0,1d0,0d0/))
      parameter (zaxis = (/0d0,0d0,1d0/))
      logical, save :: fixed_lepton_beam = .true.


C -   N.B: To check the phase space volume is correct,
C -   replace res(j) in sigborn.f by res(j)=kn_jacborn
C -   then check that the cross section is equal to,
C -   hc^2/(8*pi) = 1.54929*10^7.
c~       print*, 'xborn', xborn
      ! Initialise stuff
      if(ini) then
C -   Set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0

!     Read in bounds on x,y,Q2
         Q2min = (powheginput('Qmin'))**2
         Q2max = (powheginput('Qmax'))**2
         xmin = powheginput('xmin')
         xmax = powheginput('xmax')
         ymin = powheginput('ymin')
         ymax = powheginput('ymax')

         
!     Set centre of mass energy
         El = kn_beams(0,1)
         Eh = kn_beams(0,2)
         s = 4d0 * Eh * El

!     Sanity checks
!         if(Q2min.lt.0d0) stop 'Q2min negative'
         if(Q2min.lt.1d0) stop 'Qmin must be at least 1 GeV'
         if(xmin.lt.0d0) stop 'xmin negative'
         if(ymin.lt.0d0) stop 'ymin negative'
         if(Q2min.gt.Q2max) stop 'Q2min > Q2max'
         if(xmin.gt.xmax) stop 'xmin > xmax'
         if(ymin.gt.ymax) stop 'ymin > ymax'
         if(Q2max.gt.s) Q2max = s

         if(Q2min.eq.Q2max.and.xmin.eq.xmax.and.ymin.eq.ymax) stop
     $        'Only two variables can be constrained'

         if(s*xmax*ymax.lt.Q2min .or. s*xmin*ymin.gt.Q2max) stop
     $        'No phase space available'

         if(Q2min.eq.Q2max.and.xmin.eq.xmax) then
            ymin = Q2min/(s*xmin)
            ymax = ymin
         else
            if(xmin*(ymax*s) .lt. Q2min) xmin = Q2min/(ymax*s)
            if(xmax*(ymin*s) .gt. Q2max) xmax = Q2max/(ymin*s)
         endif
!         if(s*xmin*ymin .gt. Q2min) Q2min = s*xmin*ymin
!     if(s*xmax*ymax .gt. Q2max) Q2max = s*xmax*ymax


!     Q^2 = x1 x2 y S
!     x1 = Q^2/x2/y/S

         fixed_lepton_beam =(powheginput('#fixed_lepton_beam') .ne. 0d0)
         if(fixed_lepton_beam) then
            xminl = 1d0
            xmaxl = 1d0
         else
            xminl = Q2min/xmax/ymax/s
            if(xmin*ymin .eq. 0d0 ) then
               xmaxl = 1d0
            else
               xmaxl = min(Q2max/xmin/ymin/s,1d0)
            endif
         endif
         if(xminl > xmaxl) stop
     $        'No phase space available'

         print*, '****************************************'
         print*, 'Doing DIS with'
         print*, 'xmin, xmax', xmin, xmax
         print*, 'ymin, ymax', ymin, ymax
         print*, 'Q2min, Q2max', Q2min, Q2max
         print*, 'xminl, xmaxl', xminl, xmaxl
         print*, '****************************************'

         !read(*,*)
         ini=.false.
      endif
c~       xborn = (/ 0.1d0, 0.2d0, 0.3d0/)
c~       print*, 'xborn', xborn
!     Initial jacobian pre factor
c~       jac = 1d0/(16d0*pi)
      jac = 1d0/(32d0*pi**2)                     !!!! Phi = dX dOmega /(32pi^2)
!     Generate first momentum fraction x
      ixb = 1
      if(xmin.eq.xmax) then
         x = xmin
         jac = 1d0 * jac
!      elseif(xmin.lt.ymax) then
      elseif(xmin.lt.xmax) then
         x = xmin * exp(xborn(ixb)*log(xmax/xmin))
         jac = jac * log(xmax/xmin) * x         !!!! Include the jacobian for dX /dxborn(1)
      endif

!     Now generate the lepton x
      ixb=ixb+1
      if(fixed_lepton_beam) then
         xl = xminl
         jac = 1d0 * jac
      else
!     COMMENT -- ADD IMPORTANCE SAMPLING IF NEEDED

!     x1 = Q^2/x2/y/S
         if(ymin .eq. 0) then
            xlp = xmaxl
         else
            xlp = min(Q2max/x/s/ymin, xmaxl)
         endif
         xlm = max(Q2min/x/s/ymax, xminl)
         
!xl = xlm + (xlp-xlm) * xborn(ixb)
!jac = jac * (xlp-xlm)
         call sample_x_lepton(xborn(ixb), xlm, xlp, xl, xljac)
         jac = jac * xljac
      endif   
         
      
      kn_sborn = x * xl * s
      sqrts = sqrt(kn_sborn)
      kn_xb1 = xl               ! Lepton
      kn_xb2 = x                ! Quark
      
      kn_cmpborn = 0d0
!     Initial state lepton in CM frame
      kn_cmpborn(0,1) = 1d0
      kn_cmpborn(3,1) = 1d0
!     Initial state quark in CM frame
      kn_cmpborn(0,2) =  1d0
      kn_cmpborn(3,2) = -1d0
      
!     Now generate the lepton-quark system back to back in the CM frame
!     Through generation of y
!     Check available phase space for y/x
      ymn = ymin
      ymx = ymax

      if(ymin.lt. Q2min/(s*x*xl)) ymn = Q2min/(s*x*xl)
      if(ymax.gt. Q2max/(s*x*xl)) ymx = Q2max/(s*x*xl)


!     Now generate ydis
      ixb=ixb+1
      if(ymn.eq.ymx) then
         y = ymn
         jac = 1d0 * jac
         if(Q2min .eq. Q2max) jac = jac/(x*s*xl)   !! If I fix Q2, I need the jacobian dy/dQ2
      elseif(ymn.lt.ymx) then
         y = ymn * exp(xborn(ixb)*log(ymx/ymn))
         jac = jac * log(ymx/ymn) * y         !!!! Include the jacobian for dY /dxborn(2)     
      else
         print*, ymn,ymx,x,xl,s,Q2min,Q2max
         stop 'No phase space available for y'
      endif
! For numerical stability
      if(Q2min.eq.Q2max) then
         Q2 = Q2min
      else
         Q2 = y * x * xl * s
      endif
c~       print*, 'Rolled Q2 in Born phase space', Q2
      if(Q2.lt.Q2min) then
         print*, 'Q2 too low', Q2, Q2min
      elseif(Q2.gt.Q2max) then
         print*, 'Q2 too high', Q2, Q2max
      endif
!      print*, Q2
!     theta = acos(1-2d0*y)
      cost = 1-2d0*y
      sint = sqrt(1d0 - cost**2)
      jac = 2d0 * jac ! y = (1-cos(theta))/2   !!!! Include the jacobian for dcth /dY


      ixb=ixb+1
      phi = xborn(ixb) * 2d0 * pi
      jac = jac * 2d0 * pi

!     Outgoing lepton in CM frame
      kn_cmpborn(0,3) = 1d0
      kn_cmpborn(1,3) = Cos(phi) * sint
      kn_cmpborn(2,3) = Sin(phi) * sint
      kn_cmpborn(3,3) = cost
!     Rotate around y and z axis
c~       call mrotate(yaxis,sint,cost,kn_cmpborn(1:3,3))
c~       call mrotate(zaxis,sin(phi),cos(phi),kn_cmpborn(1:3,3))
!      call mrotate(yaxis,sin(theta),cos(theta),kn_cmpborn(1:3,3))

!     Quark back to back
      kn_cmpborn(0,4) = 1d0
      kn_cmpborn(1:3,4) = -kn_cmpborn(1:3,3)
!     Now rescale to total CM energy
      kn_cmpborn = sqrts * 0.5d0 * kn_cmpborn

!     Boost into lab frame
      beta=(xl*El-x*Eh)/(x*Eh+xl*El)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,beta,kn_cmpborn(0,1),kn_pborn(0,1))

      q = kn_pborn(:,1) - kn_pborn(:,3)
      call mlab2breit2(4,q,kn_pborn,kn_bfpborn,.False.)

      kn_jacborn = jac

c$$$      call checkmomzero(nlegborn,kn_pborn)
c$$$      if(x.lt.xmin.or.x.gt.xmax) print*, 'x out of bounds', x, xmax, xmin
c$$$      if(y.lt.ymin.or.y.gt.ymax) print*, 'y out of bounds', y, ymax, ymin
c$$$      if(Q2.lt.Q2min.or.Q2.gt.Q2max) print*, 'Q2 out of bounds', Q2, Q2max, Q2min
c$$$      print*, 'DIS variables xl, x, y, Q', xl, x, y, sqrt(Q2)
c$$$      print*, 'Incoming particles in lab frame'
c$$$      print*, kn_pborn(:,1)
c$$$      print*, kn_pborn(:,2)
c$$$      print*, kn_pborn(:,3)
c$$$      print*, kn_pborn(:,4)
c$$$      print*, ''
c$$$      print*, 'Incoming particles in CM frame'
c$$$      print*, kn_cmpborn(:,1)
c$$$      print*, kn_cmpborn(:,2)
c$$$      print*, kn_cmpborn(:,3)
c$$$      print*, kn_cmpborn(:,4)
c$$$      print*, ''
c$$$      print*, 'kn_jacoborb=', kn_jacborn
c$$$      print*, 'Event generated successfully'
c$$$  read(*,*)
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real *8 fact
      real *8, save :: q2suppr
      real *8 q2
      real *8, external :: dotp, powheginput
      if(ini) then
         ini = .false.
         q2suppr = max(powheginput("#bornsuppfact"), 0d0)
         q2suppr = q2suppr**2
      endif
      q2 = 2d0*dotp(kn_cmpborn(:,1), kn_cmpborn(:,3))
      fact = q2/(q2 + q2suppr)
      end

      subroutine regular_suppression(fact)
      implicit none
      real * 8 fact
      call rmn_suppression(fact)
      end


      subroutine global_suppression(c,fact)
      implicit none
      character * 1 c
      real * 8 fact
      fact=1d0
      end

      subroutine rmn_suppression(fact)
      implicit none    
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real *8 fact
      real *8, save :: q2suppr
      real *8 q2
      real *8, external :: dotp, powheginput
      if(ini) then
         ini = .false.
         q2suppr = max(powheginput("#bornsuppfact"), 0d0)
         q2suppr = q2suppr**2
      endif
      q2 = 2d0*dotp(kn_cmpborn(:,1), kn_cmpborn(:,3))
      fact = q2/(q2 + q2suppr)
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real * 8 muf,mur
      real * 8 q(0:3), Q2,x,y,ptlep2
      logical, save :: ini
      integer, save :: runningscales
      real * 8, external:: powheginput
      real * 8, save :: boson_mass
      data ini/.true./
      if (ini) then
         runningscales = powheginput('#runningscales')
         if(powheginput("#channel_type") .eq. 3d0) then
            boson_mass = ph_Wmass
         else
            boson_mass = ph_Zmass
         endif

         ini =.false.
         return
      endif


      if(runningscales.lt.1) then
         muf = boson_mass
         mur = muf
         mur = muf
         return
      else
         if ((flg_btildepart.eq.'r').or.(flg_btildepart.eq.'R')) then
            q = kn_preal(:,1) - kn_preal(:,3)
            ptlep2 = kn_preal(1,3)**2 + kn_preal(2,3)**2  
         else
            q = kn_pborn(:,1) - kn_pborn(:,3)
            ptlep2 = kn_pborn(1,3)**2 + kn_pborn(2,3)**2  
         endif

         Q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2
         Q2 = -Q2
         y = 1d0 - ptlep2/Q2
         x = Q2/(y*kn_sbeams)

         if(runningscales.eq.1) then ! Use Q
            muf=max(sqrt(Q2),1d0)
         elseif(runningscales.eq.2) then ! Use Q2*(1-y) = ptlep**2
            muf=max(sqrt(Q2*(1d0-y)),1d0)
         elseif(runningscales.eq.3) then ! Use Q2*(1-x)/x = maxpt**2
            muf=max(sqrt(Q2*(1d0-x)/x),1d0)
         else
            stop 'Wrong runningscales argument!'
         endif
         mur = muf
      endif

      end


   
