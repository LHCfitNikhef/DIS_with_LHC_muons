!     This analysis computes the reduced cross section for DIS in bins
!     of x as a function of Q^2 and vice versa. The reduced cross section is given by
!
!     σ_r = x Q^4 / (2 π α^2 (1 + (1-y)^2)) * dσ / dx dQ^2
!
!     ie at leading order it is the F2 structure function.

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_rwl.h'
      integer xnbins, Q2nbins
      
      real * 8 Q2binsize, Q2min, Q2max, logQ2min, logQ2max
      real * 8 xbinsize, logxmin, logxmax, sbeams
      call inihists
      call setupmulti(rwl_num_weights+1)
      
      ! Bins in Q^2 and as a function of x
      xnbins = 50
      xbinsize = 0.5d0 / dble(xnbins)
      call bookupeqbins('sigma_r_Q2_5.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_7.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_9.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_11.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_13.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_16.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_20.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_32.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_40.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_50.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_65.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_85.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_110.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_140.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_185.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_240.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_310.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_410.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_530.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_710.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_900.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_1300.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_1800.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_2500.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_3500.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_15000.0',xbinsize,0d0,0.5d0)
      
      ! Bins in Q^2 and as a function of log(x)
      logxmax = log(0.5d0)!0d0
      logxmin = -12d0
      xbinsize = (logxmax - logxmin) / dble(xnbins)
      call bookupeqbins('log_sigma_r_Q2_5.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_7.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_9.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_11.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_13.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_16.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_20.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_32.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_40.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_50.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_65.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_85.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_110.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_140.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_185.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_240.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_310.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_410.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_530.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_710.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_900.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_1300.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_1800.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_2500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_3500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_15000.0',xbinsize,logxmin,logxmax)

      Q2nbins = 50
!     Bins in x and as a function of log(Q^2)
      sbeams = 4d0 * ebmup(1) * ebmup(2)
      logQ2min = log(3d0)/log(10d0)
      logQ2max = log(sbeams)/log(10d0)
      Q2binsize = (logQ2max - logQ2min) / dble(Q2nbins)
!      print*, sbeams, logQ2min, logQ2max, Q2binsize
!      stop
      
      call bookupeqbins('log_sigma_r_x_0.00032',Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0005', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0008', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0013', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.002',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0032', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.005',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.008',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.013',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.02',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.032',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.05',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.08',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.13',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.2',    Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.32',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.5',    Q2binsize,logQ2min,logQ2max)

 
      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'nlegborn.h'
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_weights.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwl.h'
      include 'PhysPars.h'
      real * 8 dsig(weights_max)

      logical ini
      data ini/.true./
      save ini
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=128)
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer nlep, npartons                         !lepton, quark, lepton initialstate, parton initialstate 
      real * 8 plephard(4),plep(4,maxtrack),pquark(4,maxtrack),plis(4),ppis(4),q(4)
!     momenta of lepton, quark, lepton initialstate, parton is, momentum transfere
      real * 8 plep03(0:3,maxtrack),pquark03(0:3,maxtrack),plis03(0:3),ppis03(0:3)
     $     ,q03(0:3) 
      real * 8 Q2,sbeams,x,y, delx !momentum transfer variable Q2=-q^2
      real * 8, parameter :: alphaem = 1d0/137d0 ! The lhef_analysis does not fill the PhysPars.h apparently

      integer i

      real * 8, external :: phepdot,getpseudorapidity0

      real *8 sigma_r

      dsig=0d0
      call multi_plot_setup(dsig0,dsig, weights_max)
!      print*, 'weights', dsig
      
      if(ini) then         
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         ini=.false.
      endif

!      if (nhep == 4) return
      nlep = 0
      npartons = 0
      pquark = 0d0
      plep = 0d0
      ppis = 0d0
      plis = 0d0
      if(whcprg.eq.'NLO   '.or.whcprg.eq.'LHE   '.or.whcprg.eq.'PYTHIA')
     $     then
         do i=1,nhep
!     Initial states
            if(isthep(i).eq.-1.or.isthep(i).eq.21) then
               if(abs(idhep(i)).eq.11) then
                  plis(1:4) = phep(1:4,i)
               else if(abs(idhep(i)).le.5.or.idhep(i).eq.21) then
                  ppis(1:4) = phep(1:4,i)
               else
                  stop 'Wrong initial state'
               endif
!     Final states
            else if(isthep(i).eq.1) then
               if(abs(idhep(i)).eq.11) then
                  nlep = nlep + 1
                  plep(1:4,nlep) = phep(1:4,i)
                  ! For now everything not an electron is assumed a
                  ! parton. If we run with QED shower we have to be
                  ! careful
               else 
                  npartons = npartons + 1
                  pquark(1:4,npartons) = phep(1:4,i)
               endif
            endif
         enddo
      else
         print*, 'Wrong program', whcprg
         stop
      endif


      if(nlep.ne.1) then
         print*, nlep
         stop 'Wrong number of leptons'
      endif

      sbeams = 4d0 * ebmup(1) * ebmup(2)
      
      q(:) = plis(:) - plep(:,nlep)
      
      Q2 = phepdot(q,q)
      Q2 = -Q2

      if(Q2.lt.0d0 .or. Q2.gt.sbeams) then
         print*, 'Q2', Q2
         stop 'Q2 out of bounds'
      endif
      
      y = phepdot(ppis,q) / phepdot(ppis,plis)

      if(y.gt.1d0 .or. y.lt.0d0) stop 'y out of bounds'

      x = Q2 / (sbeams * y)

      if(x.gt.1d0 .or. x.lt.0d0) stop 'x out of bounds'

!     Reduced cross section
      sigma_r =  x * Q2**2
     $     / (2d0 * pi * alphaem**2 * (1d0 + (1d0 - y)**2))
      sigma_r = sigma_r / 1d6 ! convert to mb
!      print*, 'sbeams, Q2, y, x', sbeams, Q2, y, x, sigma_r
!     Now we fill histograms first in Q2 binning
      if(Q2.lt.4d0) then
         return
      elseif(Q2.lt.5d0) then
         call filld('sigma_r_Q2_5.0',    x, sigma_r * dsig / 1d0)
         call filld('log_sigma_r_Q2_5.0',    log(x), sigma_r / x * dsig /  1d0)
      elseif(Q2.lt.7d0) then
         call filld('sigma_r_Q2_7.0',    x, sigma_r * dsig / 2d0)
         call filld('log_sigma_r_Q2_7.0',    log(x), sigma_r / x * dsig /  2d0)
      elseif(Q2.lt.9d0) then
         call filld('sigma_r_Q2_9.0',    x, sigma_r * dsig / 2d0)
         call filld('log_sigma_r_Q2_9.0',    log(x), sigma_r / x * dsig /  2d0)
      elseif(Q2.lt.11d0) then
         call filld('sigma_r_Q2_11.0',    x, sigma_r * dsig / 2d0)
         call filld('log_sigma_r_Q2_11.0',    log(x), sigma_r / x * dsig /  2d0)
      elseif(Q2.lt.13d0) then
         call filld('sigma_r_Q2_13.0',    x, sigma_r * dsig / 2d0)
         call filld('log_sigma_r_Q2_13.0',    log(x), sigma_r / x * dsig /  2d0)
      elseif(Q2.lt.16d0) then
         call filld('sigma_r_Q2_16.0',    x, sigma_r * dsig / 3d0)
         call filld('log_sigma_r_Q2_16.0',    log(x), sigma_r / x * dsig /  3d0)
      elseif(Q2.lt.20d0) then
         call filld('sigma_r_Q2_20.0',    x, sigma_r * dsig / 4d0)
         call filld('log_sigma_r_Q2_20.0',    log(x), sigma_r / x * dsig /  4d0)
      elseif(Q2.lt.32d0) then
         call filld('sigma_r_Q2_32.0',    x, sigma_r * dsig / 12d0)
         call filld('log_sigma_r_Q2_32.0',    log(x), sigma_r / x * dsig /  12d0)
      elseif(Q2.lt.40d0) then
         call filld('sigma_r_Q2_40.0',    x, sigma_r * dsig / 8d0)
         call filld('log_sigma_r_Q2_40.0',    log(x), sigma_r / x * dsig /  8d0)
      elseif(Q2.lt.50d0) then
         call filld('sigma_r_Q2_50.0',    x, sigma_r * dsig / 10d0)
         call filld('log_sigma_r_Q2_50.0',    log(x), sigma_r / x * dsig /  10d0)
      elseif(Q2.lt.65d0) then
         call filld('sigma_r_Q2_65.0',    x, sigma_r * dsig / 15d0)
         call filld('log_sigma_r_Q2_65.0',    log(x), sigma_r / x * dsig /  15d0)
      elseif(Q2.lt.85d0) then
         call filld('sigma_r_Q2_85.0',    x, sigma_r * dsig / 20d0)
         call filld('log_sigma_r_Q2_85.0',    log(x), sigma_r / x * dsig /  20d0)
      elseif(Q2.lt.110d0) then
         call filld('sigma_r_Q2_110.0',    x, sigma_r * dsig / 25d0)
         call filld('log_sigma_r_Q2_110.0',    log(x), sigma_r / x * dsig /  25d0)
      elseif(Q2.lt.140d0) then
         call filld('sigma_r_Q2_140.0',    x, sigma_r * dsig / 30d0)
         call filld('log_sigma_r_Q2_140.0',    log(x), sigma_r / x * dsig /  30d0)
      elseif(Q2.lt.185d0) then
         call filld('sigma_r_Q2_185.0',    x, sigma_r * dsig / 45d0)
         call filld('log_sigma_r_Q2_185.0',    log(x), sigma_r / x * dsig /  45d0)
      elseif(Q2.lt.240d0) then
         call filld('sigma_r_Q2_240.0',    x, sigma_r * dsig / 55d0)
         call filld('log_sigma_r_Q2_240.0',    log(x), sigma_r / x * dsig /  55d0)
      elseif(Q2.lt.310d0) then
         call filld('sigma_r_Q2_310.0',    x, sigma_r * dsig / 70d0)
         call filld('log_sigma_r_Q2_310.0',    log(x), sigma_r / x * dsig /  70d0)
      elseif(Q2.lt.410d0) then
         call filld('sigma_r_Q2_410.0',    x, sigma_r * dsig / 100d0)
         call filld('log_sigma_r_Q2_410.0',    log(x), sigma_r / x * dsig /  100d0)
      elseif(Q2.lt.530d0) then
         call filld('sigma_r_Q2_530.0',    x, sigma_r * dsig / 120d0)
         call filld('log_sigma_r_Q2_530.0',    log(x), sigma_r / x * dsig /  120d0)
      elseif(Q2.lt.710d0) then
         call filld('sigma_r_Q2_710.0',    x, sigma_r * dsig / 180d0)
         call filld('log_sigma_r_Q2_710.0',    log(x), sigma_r / x * dsig /  180d0)
      elseif(Q2.lt.900d0) then
         call filld('sigma_r_Q2_900.0',    x, sigma_r * dsig / 190d0)
         call filld('log_sigma_r_Q2_900.0',    log(x), sigma_r / x * dsig /  190d0)
      elseif(Q2.lt.1300d0) then
         call filld('sigma_r_Q2_1300.0',    x, sigma_r * dsig / 400d0)
         call filld('log_sigma_r_Q2_1300.0',    log(x), sigma_r / x * dsig /  400d0)
      elseif(Q2.lt.1800d0) then
         call filld('sigma_r_Q2_1800.0',    x, sigma_r * dsig / 500d0)
         call filld('log_sigma_r_Q2_1800.0',    log(x), sigma_r / x * dsig /  500d0)
      elseif(Q2.lt.2500d0) then
         call filld('sigma_r_Q2_2500.0',    x, sigma_r * dsig / 700d0)
         call filld('log_sigma_r_Q2_2500.0',    log(x), sigma_r / x * dsig /  700d0)
      elseif(Q2.lt.3500d0) then
         call filld('sigma_r_Q2_3500.0',    x, sigma_r * dsig / 1000d0)
         call filld('log_sigma_r_Q2_3500.0',    log(x), sigma_r / x * dsig /  1000d0)
      elseif(Q2.lt.15000d0) then
         call filld('sigma_r_Q2_15000.0',    x, sigma_r * dsig / 11500d0)
         call filld('log_sigma_r_Q2_15000.0',    log(x), sigma_r / x * dsig /  11500d0)
      endif

      if(x.lt.0.0002d0) then
         return
      elseif(x.lt.0.00032d0) then
         delx = 0.00032d0 - 0.0002d0
         call filld('log_sigma_r_x_0.00032',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.0005d0) then
         delx = 0.0005d0 - 0.00032d0
         call filld('log_sigma_r_x_0.0005',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.0008d0) then
         delx = 0.0008d0 - 0.0005d0
         call filld('log_sigma_r_x_0.0008',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.0013d0) then
         delx = 0.0013d0 - 0.0008d0 
         call filld('log_sigma_r_x_0.0013',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.002d0) then
         delx = 0.0020d0 - 0.0013d0
         call filld('log_sigma_r_x_0.002',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.0032d0) then
         delx = 0.0032d0 - 0.0020d0
         call filld('log_sigma_r_x_0.0032',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.005d0) then
         delx = 0.005d0 - 0.0032d0
         call filld('log_sigma_r_x_0.005',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.008d0) then
         delx = 0.008d0 - 0.005d0
         call filld('log_sigma_r_x_0.008',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.013d0) then
         delx = 0.013d0 - 0.008d0
         call filld('log_sigma_r_x_0.013',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.02d0) then 
         delx = 0.02d0 - 0.013d0
        call filld('log_sigma_r_x_0.02',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.032d0) then
         delx = 0.032d0 - 0.02d0
         call filld('log_sigma_r_x_0.032',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.05d0) then
         delx = 0.05d0 - 0.032d0
         call filld('log_sigma_r_x_0.05',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.08d0) then
         delx = 0.08d0 - 0.05d0
         call filld('log_sigma_r_x_0.08',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.13d0) then
         delx = 0.13d0 - 0.08d0
         call filld('log_sigma_r_x_0.13',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.2d0) then
         delx = 0.2d0 - 0.13d0
         call filld('log_sigma_r_x_0.2',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.32d0) then
         delx = 0.32d0 - 0.2d0
         call filld('log_sigma_r_x_0.32',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      elseif(x.lt.0.5d0) then
         delx = 0.5d0 - 0.32d0
         call filld('log_sigma_r_x_0.5',log(Q2)/log(10d0), sigma_r * dsig / Q2 / delx)
      endif


      end



      function phepDot(p_A,p_B)
      implicit none
      real * 8  phepDot
      real * 8  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end

      function kt2(p)
      implicit none
      real * 8 kt2, p(0:3)

      kt2 = p(1)**2 + p(2)**2
      end

      function eta(p)
      implicit none
      real * 8 eta, p(0:3), normp, norm

      normp = norm(p)
      if(normp.gt.p(3)) then
         eta = 0.5d0 * log((normp + p(3)) / (normp - p(3)))
      else
         eta = sign(1d100,p(3)) 
      endif
      end

      function norm(p)
      implicit none
      real * 8 norm, p(0:3)

      norm = p(1)**2 + p(2)**2 + p(3)**2
      norm = sqrt(max(0d0,norm))
      end

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end
      
      function getpseudorapidity0(p)
      implicit none
      real * 8 p(0:3),getpseudorapidity0,norm      
      getpseudorapidity0=0.5d0*log((norm(p)+p(3))/(norm(p)-p(3)))
      end

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function azi(p)
      implicit none
      real * 8 pi
      parameter(pi = 3.141592653589793D0)
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end
      

