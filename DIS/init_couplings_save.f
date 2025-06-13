      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include "pwhg_physpar.h"
      include 'vtype.h'
      include 'pwhg_flg.h'

      real * 8 masswindow
      real * 8 powheginput
      external powheginput

      integer i,j

      integer vtype

      logical verbose
      parameter(verbose=.true.)

      
      real *8 q2cutB, q2cutR, q2cutA
      common /q2cut/q2cutB, q2cutR, q2cutA



      q2cutB = max(1d0, powheginput("#q2cutB"))
      q2cutR = max(1d0, powheginput("#q2cutR"))

      flg_tiny_alphas = (powheginput("#tiny_alphas")==1d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c decide which types of contributions are considered for NC:
c 1: photon exchange only / 2: Z exchange only / 3: photon+Z exchange
      vtype =  powheginput("vtype")

      aonly = .false.
      zonly = .false.
      if (vtype.eq.1) then
         aonly = .true.
      elseif (vtype.eq.2) then
         zonly = .true.
      endif !vtype   

      ph_Hmass = powheginput("#hmass")
      if(ph_Hmass.lt.0d0) ph_Hmass  = 125d0
      ph_Hwidth = powheginput("#hwidth")
      if (ph_Hwidth.lt.0d0) ph_Hwidth = 4.14d-03
      ph_topmass = powheginput("#topmass")
      if(ph_topmass<0d0) ph_topmass = 172.5d0
      ph_Zmass = powheginput("#Zmass")
      if(ph_Zmass<0d0) ph_Zmass  = 91.1876d0
      ph_Zwidth = powheginput("#Zwidth")
      if(ph_Zwidth<0d0) ph_Zwidth =  2.4952d0
      ph_Wmass = powheginput("#Wmass")
      if(ph_Wmass<0d0) ph_Wmass  = 80.398d0
      ph_Wwidth = powheginput("#Wwidth")
      if(ph_Wwidth<0d0) ph_Wwidth = 2.141d0
      ph_alphaem = powheginput("#alphaem")
      if(ph_alphaem<00d0) ph_alphaem = 1d0/137d0
      ph_sthw2 = 1d0 - (ph_Wmass/ph_Zmass)**2 

c     number of light flavors
      st_nlight = 5

      ph_CKM(1,1)=0.9748 	
      ph_CKM(1,2)=0.2225  	 
      ph_CKM(1,3)=0.0036  	
      ph_CKM(2,1)=0.2225  	
      ph_CKM(2,2)=0.9740 	
      ph_CKM(2,3)=0.041	
      ph_CKM(3,1)=0.009    
      ph_CKM(3,2)=0.0405   
      ph_CKM(3,3)=0.9992

c     initialize CKM with flavor indexes
      call inizialize_ph_CKM_matrix

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2
      ph_Wmass2 = ph_Wmass**2
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_WmWw = ph_Wmass * ph_Wwidth
      masswindow = 30
      ph_Hmass2low=max(0d0,ph_Hmass-masswindow*ph_Hwidth)
      ph_Hmass2low=ph_Hmass2low**2
      ph_Hmass2high=(ph_Hmass+masswindow*ph_Hwidth)**2
      
      ph_unit_e = sqrt(4*pi*ph_alphaem)



 
c     Set here lepton and quark masses for momentum reshuffle in the LHE event file
      physpar_ml(1) = 0.51099891d-3
      physpar_ml(2) = 0.1056583668d0
      physpar_ml(3) = 1.77684d0
      physpar_mq(1) = 0.33d0    ! down
      physpar_mq(2) = 0.33d0    ! up
      physpar_mq(3) = 0.50d0    ! strange
      physpar_mq(4) = 1.50d0    ! charm
      physpar_mq(5) = 4.5d0     ! bottom
      
         

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) 'H mass = ',ph_Hmass
      write(*,*) 'H width = ',ph_Hwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
      write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
      write(*,*) 'CKM matrix' 
      do i=1,3
         write(*,*) (ph_CKM(i,j),j=1,3)
      enddo
      write(*,*) '*************************************'
      endif      


c convert couplings into format needed by EW matrixelements:
      call coup_powheg_to_vbfnlo

      end


      subroutine inizialize_ph_CKM_matrix
      implicit none     
      include 'PhysPars.h'  
      integer i,j
      do i=1,6
         do j=1,6
            ph_CKM_matrix(i,j) = 0
         enddo
      enddo
      ph_CKM_matrix(1,2) = ph_CKM(1,1)
      ph_CKM_matrix(1,4) = ph_CKM(2,1)
      ph_CKM_matrix(1,6) = ph_CKM(3,1)
      ph_CKM_matrix(2,3) = ph_CKM(1,2)
      ph_CKM_matrix(2,5) = ph_CKM(1,3)
      ph_CKM_matrix(3,4) = ph_CKM(2,2)
      ph_CKM_matrix(3,6) = ph_CKM(3,2)
      ph_CKM_matrix(4,5) = ph_CKM(2,3)
      ph_CKM_matrix(5,6) = ph_CKM(3,3)
      do i=1,6
         do j=i+1,6
            ph_CKM_matrix(j,i) = ph_CKM_matrix(i,j)
         enddo
      enddo
      end
