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


      flg_polarised_dis = (powheginput("#polarisation") .eq. 1d0)
     
      if(flg_polarised_dis) then
         i = abs(powheginput("ih1"))
         if(i .ne. 11 .and. i .ne. 13 .and. i .ne. 15) then
            write(*,*) "Polarised DIS available only for charged ",
     $           "lepton beams. Aborting "
            stop
         endif
         write(*,*) "Doing polarised DIS: (++) + (--) - (+-) - (-+)"
      else
         write(*,*) "Doing upolarised DIS"
      endif

      flg_tiny_alphas = (powheginput("#tiny_alphas")==1d0)

      q2cutB = 1d0 
      q2cutR = 1d0 
      
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
      ph_Hwidth = powheginput("#hwidth")
      
      ph_topmass = 172.5d0

c     if one of two parameters is missing, use the default ones
      if ((ph_Hmass.lt.0d0).or.(ph_Hwidth.lt.0d0)) then
         ph_Hmass  = 125d0
         ph_Hwidth = 4.14d-03
      endif

      write(*,*) '**************************************'
      write(*,*) '**************************************'
      write(*,*) 'Higgs boson mass  = ',ph_Hmass
      write(*,*) 'Higgs boson width = ',ph_Hwidth
      write(*,*) '**************************************'
      write(*,*) '**************************************'

      !Allow to over write the hardcoded values
      ! Can in principle provide the pole masses (rather than OS) in powheg.input
      ! m_pole = m_os / sqrt(1 + (G_os/M_os)^2), G_pole = G_os / sqrt(1 + (G_os/M_os)^2)
      ph_Zmass = powheginput("#Zmass")
      if(ph_Zmass<0d0) ph_Zmass  = 91.1876d0
      ph_Zwidth = powheginput("#Zwidth")
      if(ph_Zwidth<0d0) ph_Zwidth =  2.4952d0
      ph_Wmass = powheginput("#Wmass")
      if(ph_Wmass<0d0) ph_Wmass  = 80.398d0
      ph_Wwidth = powheginput("#Wwidth")
      if(ph_Wwidth<0d0) ph_Wwidth = 2.141d0
      ! For the CC process, it is preferable to include top mass and 
      ! fermion mass effects. So use:
      ! alpha_GF = sqrt(2) / M_PI * gf * fabs( MW2C * SW2_OS )
      ph_alphaem = powheginput("#alphaem")
      if(ph_alphaem<0d0) ph_alphaem = 1d0/137d0
      ! In principle use the all-order on-shell relation for S2W
      ! At LO in EW, for better accuracy one can instead fix s2w or 
      ! include rho-parameter corrections
      ph_sthw2 = powheginput("#s2w")
      if(ph_sthw2<0d0) ph_sthw2 = 1d0 - (ph_Wmass/ph_Zmass)**2      
      
c     number of light flavors
      st_nlight = 5

      if(powheginput("#third_gen_diag_CKM").eq.1d0) then
         ph_CKM(3,:) = 0d0
         ph_CKM(:,3) = 0d0
         ph_CKM(1,2) = 0.224561d0
         ph_CKM(1,1) = 0.97446d0
         ph_CKM(2,1) = ph_CKM(1,2)
         ph_CKM(2,2) = ph_CKM(1,1)
      else
         ph_CKM(1,1) = powheginput('#CKM_Vud')
         ph_CKM(1,2) = powheginput('#CKM_Vus')
         ph_CKM(1,3) = powheginput('#CKM_Vub')
         ph_CKM(2,1) = powheginput('#CKM_Vcd')
         ph_CKM(2,2) = powheginput('#CKM_Vcs')
         ph_CKM(2,3) = powheginput('#CKM_Vcb')
         ph_CKM(3,1) = powheginput('#CKM_Vtd')
         ph_CKM(3,2) = powheginput('#CKM_Vts')
         ph_CKM(3,3) = powheginput('#CKM_Vtb')

         
         if(ph_CKM(1,1)<0) ph_CKM(1,1)=0.9748 	
         if(ph_CKM(1,1)<0) ph_CKM(1,2)=0.2225  	 
         if(ph_CKM(1,1)<0) ph_CKM(1,3)=0.0036  	
         if(ph_CKM(1,1)<0) ph_CKM(2,1)=0.2225  	
         if(ph_CKM(1,1)<0) ph_CKM(2,2)=0.9740 	
         if(ph_CKM(1,1)<0) ph_CKM(2,3)=0.041	
         if(ph_CKM(1,1)<0) ph_CKM(3,1)=0.009    
         if(ph_CKM(1,1)<0) ph_CKM(3,2)=0.0405   
         if(ph_CKM(1,1)<0) ph_CKM(3,3)=0.9992
      endif


c     initialize CKM with flavor indexes
      call inizialize_ph_CKM_matrix

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2
      ph_Wmass2 = ph_Wmass**2
      ph_Hmass2 = ph_Hmass**2
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_HmHw = ph_Hmass * ph_Hwidth


c     set mass windows around the resonance peak 
c     It is used in the generation of the Born phase space
      masswindow = 30
c      ph_Zmass2low=(ph_Zmass-masswindow*ph_Zwidth)**2
c      ph_Zmass2high=(ph_Zmass+masswindow*ph_Zwidth)**2
      ph_Hmass2low=max(0d0,ph_Hmass-masswindow*ph_Hwidth)
      ph_Hmass2low=ph_Hmass2low**2
      ph_Hmass2high=(ph_Hmass+masswindow*ph_Hwidth)**2
c      ph_Hmass2low=0d0
c      ph_Hmass2high=kn_sbeams/4
    
      ph_unit_e = sqrt(4*pi*ph_alphaem)



      if(powheginput("#masslesslhe") == 1d0) then
         physpar_ml = 0d0
         physpar_mq = 0d0
      else    
c     Set here lepton and quark masses for momentum reshuffle in the LHE event file
         physpar_ml(1) = powheginput("#electron_mass")
         physpar_ml(2) = powheginput("#muon_mass")
         physpar_ml(3) = powheginput("# tauon_mass")
         if(physpar_ml(1) .lt. 0d0) physpar_ml(1) = 0.51099891d-3
         if(physpar_ml(2) .lt. 0d0) physpar_ml(2) = 0.1056583668d0
         if(physpar_ml(3) .lt. 0d0) physpar_ml(3) = 1.77684d0
         
         physpar_mq(1) = 0.33d0 ! down
         physpar_mq(2) = 0.33d0 ! up
         physpar_mq(3) = 0.50d0 ! strange
         physpar_mq(4) = powheginput("#charm_mass")
         physpar_mq(5) = powheginput("#bottom_mass")
         if(physpar_mq(4) .lt. 0d0) physpar_mq(4) = 1.50d0 ! charm
         if(physpar_mq(5) .lt. 0d0) physpar_mq(5) = 4.5d0 ! bottom
      endif
         

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
