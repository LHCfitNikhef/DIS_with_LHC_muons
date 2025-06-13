!  The next subroutines, open some histograms and prepare them 
!      to receive data 
!  You can substitute these  with your favourite ones
!  init   :  opens the histograms
!  topout :  closes them
!  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_rwl.h'

      integer i,n
      parameter (n = 8)
      real *8, parameter ::  Qbins(1:n) =
     $     (/14d0, 16d0, 20d0, 30d0, 50d0, 70d0, 100d0, 200d0/)
      character(len=8), parameter ::  Qchar(1:n) =
     $     (/'____0-14','___14-16','___16-20','___20-30',
     $     '___30-50','___50-70','__70-100', '_100-200'/) 

      call inihists
      call setupmulti(rwl_num_weights+1)

      call bookupeqbins('sigtot',1d0,0d0,1d0)
      call bookupeqbins('Q2', (10000d0/20d0), 0d0, 10000d0)
      call bookupeqbins('y',  0.05d0 , 0d0, 1d0)
      call bookupeqbins('x',  0.05d0 , 0d0, 1d0)

      call bookupeqbins('logbroad',  (12d0/64d0), -10d0, 2d0)
      call bookupeqbins('brd_zE',0.2d0,0d0,1d0)
      call bookupeqbins('brd_zE_log',1d0,-5d0,0d0)
      call bookupeqbins('tau',0.2d0,0d0,1d0)
      call bookupeqbins('tau_log',1d0,-5d0,0d0)
      call bookupeqbins('tau_C',0.2d0,0d0,1d0)
      call bookupeqbins('tau_C_log',1d0,-5d0,0d0)
      call bookupeqbins('rho',0.2d0,0d0,1d0)
      call bookupeqbins('rho_log',1d0,-5d0,0d0)
      call bookupeqbins('cprE',0.2d0,0d0,1d0)
      call bookupeqbins('cprE_log',1d0,-5d0,0d0)

      
      do i = 2,n
         call bookupeqbins('sig_after_cuts'//trim(Qchar(i)),1d0,0d0,1d0)
         call bookupeqbins('tau'//trim(Qchar(i)),0.05d0,0d0,1d0)
         call bookupeqbins("tau_C"//trim(Qchar(i)),0.05d0,0d0,0.5d0)      
         call bookupeqbins("brd_zE"//trim(Qchar(i)),0.05d0,0d0,0.5d0)
         call bookupeqbins("rho"//trim(Qchar(i)),0.03125d0,0d0,0.25d0)
         call bookupeqbins("cprE"//trim(Qchar(i)),0.03125d0,0d0,0.25d0)
      end do

      end
     
      subroutine analysis(dsig0)
      use observable_consts
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
      real * 8 dsig(weights_max)

      logical ini
      data ini/.true./
      save ini
!     binsize
!     we need to tell to this analysis file which program is running it
      integer   maxtrack,maxjet
      integer i,n
      parameter (n = 8)
      parameter (maxtrack=2048,maxjet=128)
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer nlep, npartons                         !lepton, quark, lepton initialstate, parton initialstate 
      real * 8 plep(4,maxtrack), plephard(4),pquark(4,maxtrack),plis(4),ppis(4),q(4)              !momenta of lepton, quark, lepton initialstate, parton is, momentum transfere
      real * 8 q03(0:3), q03breit(0:3), El_hardest
      real * 8 pbreit(0:3,maxtrack),plab(0:3,maxtrack)
      real * 8 pevsh(4,maxtrack)
      real * 8 y,x,Q2,sbeams               !xs. inv mass of incoming parton-electron, momentum transfer variable Q2=-q^2
      real * 8 broad
      real * 8 brd_zE,tau, tau_c, rho, cprE, tau_q
      real * 8 qevsh(4,0:4), obsf77, kt2, norm, Ecur, pcur
      external obsf77
      real *8, parameter ::  Qbins(1:n) =
     $     (/14d0, 16d0, 20d0, 30d0, 50d0, 70d0, 100d0, 200d0/)
      character(len=8), parameter ::  Qchar(1:n) =
     $     (/'____0-14','___14-16','___16-20','___20-30',
     $     '___30-50','___50-70','__70-100', '_100-200'/) 

      real * 8, external :: phepdot


      
      dsig=0d0
      call multi_plot_setup(dsig0,dsig, weights_max)
      if(sum(abs(dsig)).eq.0) return

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

!      if (rwl_type <= 1) return

      if(nlep<1) stop 'Wrong number of leptons'

      El_hardest = 0d0
      do i = 1, nlep
C         write(*,*) 'plep(:,i)',plep(:,i)
         if (plep(4,i) > El_hardest ) then
            El_hardest = plep(4,i)
            plephard(1:4) = plep(1:4,i) 
         endif
      end do

      sbeams = 4d0 * ebmup(1) * ebmup(2)

      q(:) = plis(:) - plephard(:)
      
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
      if ((Q2<196d0 .or. Q2>40000) .or. (y<0.1d0 .or. y>0.7d0)) return
      
      call filld('sigtot',0.5d0,dsig)
      call filld('Q2', Q2, dsig)
      call filld('x', x, dsig)
      call filld('y', y, dsig)

      do i=2,8 
         if (Qbins(i-1)< sqrt(Q2) .and. sqrt(Q2)< Qbins(i)) then
            call filld('sig_after_cuts'//trim(Qchar(i)),0.5d0,dsig)
         endif
      enddo

      q03(:) = cshift(q,-1) 
      pbreit = 0d0
      plab = 0d0
      plab(:,1) = cshift(plis,-1)
      plab(:,2) = cshift(ppis,-1)
      plab(:,3:2+nlep) = cshift(plep(:,1:nlep),-1,dim=1)
      ! do i=1,npartons
      plab(:,3+nlep:2+nlep+npartons) = cshift(pquark(:,1:npartons),-1,dim=1)
      ! enddo
      call mlab2breit2(2+nlep+npartons,q03,plab,pbreit,.false.)
      call mlab2breit2(1,q03,q03,q03breit,.false.)

      pevsh(:,1:npartons) = cshift(pbreit(:,3+nlep:2+nlep+npartons),1, DIM=1)

      qevsh = 0d0
      qevsh(:,0) = cshift(q03breit,1)


      broad = 0d0
      Ecur = 0d0
      pcur = 0d0
      do i =1,npartons
         if (sign(1d0,pbreit(3,2+nlep+i)) == sign(1d0,qevsh(3,0))) then
            broad = broad + sqrt(kt2(pbreit(:,2+nlep+i)))
            Ecur = Ecur + pbreit(0,2+nlep+i)
            pcur = pcur + norm2(pbreit(1:3,2+nlep+i))
         endif
      enddo
      if (Ecur<dsqrt(Q2)/10d0) return
      broad = broad / (2d0*pcur)
!     -- Thrust with respect to the photon axis (normalized to E in current hemishepere)
      tau    = obsf77(iobs_tau_zE,npartons,pevsh,qevsh)
!     -- Thrust with respect to the thrust axis (normalized to E)  
      tau_c  = obsf77(iobs_tau_tE,npartons,pevsh,qevsh)
!     -- Broadening with respect to the photon axis (normalized to E)       
      brd_zE = obsf77(iobs_brd_zE,npartons,pevsh,qevsh)

!     -- hemisphere mass (normalized to E in current hemishepere)       
      rho    = obsf77(iobs_rho_E,npartons,pevsh,qevsh)

!     -- C-parameter in current hemisphere (normalized to E in current hemishepere)             
      cprE    = obsf77(iobs_cpr_E,npartons,pevsh,qevsh)

!     -- Thrust with respect to the photon axis (normalized to Q/2)
      tau_q  = obsf77(iobs_tau_zQ,npartons,pevsh,qevsh)      

      call filld('brd_zE', brd_zE, dsig)
      call filld('brd_zE_log', log(brd_zE), dsig)
      call filld('tau', tau, dsig)
      call filld('tau_log', log(tau), dsig)
      call filld('tau_C', tau_C, dsig)
      call filld('tau_C_log', log(tau_C), dsig)
      call filld('rho', rho, dsig)
      call filld('rho_log', log(rho), dsig)
      call filld('cprE', cprE, dsig)
      call filld('cprE_log', log(cprE), dsig)

      call filld('logbroad',log(broad),dsig)


      
      do i=2,8 
         if (Qbins(i-1)< sqrt(Q2) .and. sqrt(Q2)< Qbins(i)) then
            call filld('tau'//trim(Qchar(i)),tau,dsig)
            call filld('tau_C'//trim(Qchar(i)),tau_C,dsig)
            call filld('brd_zE'//trim(Qchar(i)),brd_Ze,dsig) 
            call filld('rho'//trim(Qchar(i)),rho,dsig)
            call filld('cprE'//trim(Qchar(i)),cprE,dsig)
         endif
      enddo
      

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
      
      ! Builds pp like jets. Should take as argument all final state
      ! partons/hadrons. Returns the jets pt ordered. 
      subroutine buildjets(n,pin,pj,njets,ptj,yj,phij)
      implicit none
      integer n
      double precision pin(0:3,n) 
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=128)

!     Output
      double precision pj(0:3,maxjet)
!     Internal
      integer mu, njets, njets_all, ntracks, ijet, j
      integer jetvec(maxtrack)
      double precision pjet(4,maxjet), pjet_all(4,maxjet)
      double precision ptrack(4,maxtrack)
      double precision ptj(maxjet),yj(maxjet),phij(maxjet)
      double precision ptj_all(maxjet),yj_all(maxjet),phi_all(maxjet)
      double precision R, ptmin_fastkt, palg
      double precision azi
      external azi


      ptrack = 0d0
      jetvec = 0
      pjet = 0d0
      pjet_all = 0d0
      njets=0
      njets_all = 0
      ntracks = n
      
      ptrack(4,1:n)=pin(0,1:n)
      ptrack(1:3,1:n)=pin(1:3,1:n)
      
      R = 0.4d0
      ptmin_fastkt = 0d0
      palg = -1d0
!      palg = 0d0
!      palg = -1d0
c     -1 is anti_kt, +1 is kt

      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin_fastkt,pjet,njets,jetvec)

      if(njets.gt.maxjet) then
         print*, 'njets out of bounds!!', njets, maxjet
         stop
      endif
      
      j=0
      pjet_all=pjet
      njets_all=njets
      pjet=0d0
      njets=0d0
      ptj = 0d0
      yj = 0d0
      
      do ijet=1,njets_all
         ptj_all(ijet) = sqrt(pjet_all(1,ijet)**2 + pjet_all(2,ijet)**2)
         call getrapidity(pjet_all(:,ijet),yj_all(ijet))
         phi_all(ijet) = azi(pjet_all(:,ijet))

!         if(ptj_all(ijet).gt.ptalljetmin.and.
!     $        abs(yj_all(ijet)).lt.yjetmax) then
            j=j+1
            pjet(:,j)=pjet_all(:,ijet)
            ptj(j) = ptj_all(ijet)
            yj(j) = yj_all(ijet)
            phij(j) = phi_all(ijet)
!         endif
      enddo
      njets=j

      pj = 0d0

      do ijet=1,njets
         do mu=1,3
            pj(mu,ijet)=pjet(mu,ijet)
         enddo
         pj(0,ijet)=pjet(4,ijet)
      enddo
            
      end

