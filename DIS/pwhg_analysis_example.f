!  The next subroutines, open some histograms and prepare them 
!      to receive data 
!  You can substitute these  with your favourite ones
!  init   :  opens the histograms
!  topout :  closes them
!  pwhgfill  :  fills the histograms with data

      subroutine init_hist
         use hist2d, only: bookupeqbpdlogdd, init_dists
         implicit none
         include  'LesHouches.h'
         include 'pwhg_math.h'
         include 'nlegborn.h'
         include 'pwhg_kn.h'

         call init_dists
         call inihists
         
         call bookupeqbins('sigtot',1d0,0d0,1d0)
         call bookupeqbins('Q',  0.5d0, 2d0, 20d0)

         call bookupeqbins('Q2', (1000d0-25d0)/20d0, 25d0, 1000d0)
         call bookupeqbins('y',  (0.95d0-0.04d0)/20d0 , 0.04d0, 0.95d0)
         call bookupeqbins('x',  0.05d0 , 0d0, 1d0)
         call bookupeqbins('xl',  0.05d0 , 0d0, 1d0)   
         call bookupeqbins('ptl', 2d0, 0d0, 100d0)
         call bookupeqbins('etal',0.1d0,0d0,5d0)
         call bookupeqbins('ptj_overQ_breit', 0.01d0, 0d0, 0.5d0)

         call bookupeqbpdlogdd("xQ2",-3,0,10,0,4,10,110)
         call bookupeqbpdlogdd("charmed_xQ2",-3,0,10,0,4,10,110)

      end
       
      subroutine analysis(dsig0)
         use hist2d, only: filldd
         implicit none
         real * 8 dsig0
         include 'nlegborn.h'
         include 'hepevt.h'
         include 'pwhg_math.h' 
         include 'LesHouches.h'
         include 'pwhg_kn.h'
         include 'pwhg_weights.h'
         ! include 'pwhg_rad.h'
         ! include 'pwhg_rwl.h'
         logical ini
         data ini/.true./
         save ini
         real * 8 dsig(weights_max)
   
         integer   maxtrack,maxjet
         integer i,njets
         parameter (maxtrack=2048,maxjet=128)
         ! we need to tell to this analysis file which program is running it
         character * 6 WHCPRG
         common/cWHCPRG/WHCPRG
         ! data WHCPRG/'NLO   '/
         integer nlep, npartons,nunident                         !lepton, quark, lepton initialstate, parton initialstate 
         real * 8 plep(4,maxtrack), plephard(4),pquark(4,maxtrack),plis(4),ppis(4),q(4)              !momenta of lepton, quark, lepton initialstate, parton is, momentum transfere
         real * 8 El_hardest
         real * 8 plab(0:3,maxtrack), pjets(0:3,maxjet)
         real * 8 y,x,Q2, xl               !xs. inv mass of incoming parton-electron, momentum transfer variable Q2=-q^2
         real * 8 ptj(maxjet),yj(maxjet),phij(maxjet), eta1, ptjetmin, absEtaMax
         real * 8 ptl, etal
         real * 8, external :: phepdot, eta, kt2
         real * 8 refin(1:4), refout(1:4)
         real * 8 alpha, beta, ptbreit
         real * 8, save :: Eproton, Elepton, sbeams
         logical, save :: fixed_target
         real *8, external :: powheginput
         integer :: nchargedparticles,ncharmed
         real(kind=8) ::  W2

         
         interface
            logical function charged(pid)
               integer :: pid
            end function charged
            logical function is_charmed(pid)
               integer :: pid
            end function is_charmed
            subroutine getdphi(p1,p2,dphi,phi1,phi2)
            real(kind=8),intent(in)   :: p1(4),p2(4)
            real(kind=8),intent(out)  :: dphi
            real(kind=8),intent(out),optional :: phi1,phi2
            end subroutine getdphi
         end interface
       
        
         dsig=0d0
         call multi_plot_setup(dsig0,dsig, weights_max)
   
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

            Elepton = powheginput("ebeam1")
            Eproton = powheginput("ebeam2")
            fixed_target = (powheginput("#fixed_target") .eq. 1d0)

            if(fixed_target) then
               sbeams = 2*Elepton*Eproton + Eproton**2
            else
               sbeams = 4 * Elepton * Eproton
            endif


            if(whcprg.eq.'NLO') then
               if( abs( 4d0 * ebmup(1) * ebmup(2) - sbeams) .gt. 1d-4)
     $              then
                  write(*,*) "inconsistency in the calculation of ",
     $                 "the partonic center-of-mass energy in the ",
     $                 "analysis, aborting ..."
                  stop
               endif
            endif
            
         endif
   
         nlep = 0
         npartons = 0
         pquark = 0d0
         plep = 0d0
         ppis = 0d0
         plis = 0d0
         pquark = 0d0
         nunident = 0
         nchargedparticles=0
         ncharmed=0
         if(whcprg.eq.'NLO   '.or.whcprg.eq.'LHE   '
     $        .or.whcprg.eq.'PYTHIA') then
            do i=1,nhep
               !if(idhep(i).eq.11) print*, phep(1:4,i), eta(phep(:,i)), sqrt(kt2(phep(:,i)))
   !     Initial states
               if(isthep(i).eq.-1.or.isthep(i).eq.21) then
                  if(abs(idhep(i)).le.16 .and.abs(idhep(i)).ge.11 ) then
                     plis(1:4) = phep(1:4,i)
                  else if(abs(idhep(i)).le.5.or.idhep(i).eq.21) then
                     ppis(1:4) = phep(1:4,i)
                  else
                     stop 'Wrong initial state'
                  endif
   !     Final states
               else if(isthep(i).eq.1) then
C Count the charged tracks with an energy over 1 GeV.                  
                  if(charged(idhep(i)).and.phep(4,i).ge.1d0)then
                     nchargedparticles=nchargedparticles+1
                  end if
C Count the charmed tracks with an energy over 1 GeV.
                  if(is_charmed(idhep(i)).and.phep(4,i).ge.1d0)then
                     ncharmed=ncharmed+1
                  end if
C Extract the final state quarks and leptons.
                  if(abs(idhep(i)).ge.11 .and. abs(idhep(i)).le.16) then
                     nlep = nlep + 1
                     plep(1:4,nlep) = phep(1:4,i)
                  elseif (abs(idhep(i)) <= 9 .or. abs(idhep(i)) == 21 .or. abs(idhep(i)) > 100 ) then 
                     npartons = npartons + 1
                     pquark(1:4,npartons) = phep(1:4,i)
                  else
                     nunident = nunident + 1
                     ! print*, 'idhep(i)', idhep(i)
                  endif
               endif
            enddo
         else
            print*, 'Wrong program', whcprg
            stop
         endif

         if (nunident>0) print*, 'nunident', nunident
   
         if(nlep<1) stop 'Wrong number of leptons'
   
C Construct DIS kinematic quantities.

         El_hardest = 0d0
         do i = 1, nlep
            if (plep(4,i) > El_hardest ) then
               El_hardest = plep(4,i)
               plephard(1:4) = plep(1:4,i) 
            endif
         end do
   
         q(:) = plis(:) - plephard(:)
         Q2 = phepdot(q,q)
         Q2 = -Q2
         
         x = Q2 / (sbeams * y * xl)
         xl = plis(4)/Elepton

         y = phepdot(ppis,q) / phepdot(ppis,plis)      

         W2=Eproton**2+Q2/x*(1d0-x)
    
         
         call filld('sigtot',0.5d0,dsig)   

C Modify the cuts here.
         if(El_hardest.ge.100d0
     1      .and.nchargedparticles.ge.3
     2      .and.W2.ge.4d0
     3      .and.Q2.ge.1.65d0**2)then
            call filld('Q2', Q2, dsig)
            call filld('Q', sqrt(Q2), dsig)
            call filld('x', x, dsig)
            call filld('y', y, dsig)
            call filld('xl', xl, dsig)

            call filldd('xQ2',dsig,x,Q2)
            if(ncharmed.ge.1)then
               call filldd('charmed_xQ2',dsig,x,Q2)
            end if
         endif

         ptl  = sqrt(kt2(plephard(:)))
         call getrapidity(plephard(:),etal)
         call filld('ptl', ptl, dsig)
         call filld('etal', etal, dsig)


         if(npartons == 2) then
            refin  = x *ebmup(2) * (/ 0d0, 0d0, -1d0, 1d0/)
            if(fixed_target) refin = refin/2d0 ! The fake massless beam has E = mp/2d0, hence the factor 2
            if(plis(3) < 0) refin(3)=-refin(3) ! Flip the sign of z if necessary
            refout = q(:) + refin
            
            alpha = 2d0 * phepdot(pquark(:, 1), refout)/Q2
            beta  = 2d0 * phepdot(pquark(:, 1), refin)/Q2
            
            ptbreit = sqrt(alpha*beta)
         else
            ptbreit = 0d0
         endif
         call filld('ptj_overQ_breit', ptbreit, dsig)
  
      end
  
  
  
      function phepDot(p_A,p_B)
         implicit none
         real * 8  phepDot
         real * 8  p_A(4),p_B(4)
         phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)-p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end
  
      function kt2(p)
         implicit none
         real * 8 kt2, p(1:4)
   
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
      subroutine buildjets(n,pin,pj,njets,ptj,yj,phij, ptjetmin, absEtaMax)
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
         double precision ptj(maxjet),yj(maxjet),phij(maxjet), ptjetmin, absEtaMax
         double precision ptj_all(maxjet),yj_all(maxjet),phi_all(maxjet)
         double precision R, ptmin_fastkt, palg
         double precision, external :: azi, eta
   
   
         ptrack = 0d0
         jetvec = 0
         pjet = 0d0
         pjet_all = 0d0
         njets=0
         njets_all = 0
         ntracks = n
         
         ptrack(4,1:n)=pin(0,1:n)
         ptrack(1:3,1:n)=pin(1:3,1:n)
         
         R = 0.8d0
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
  
          if(ptj_all(ijet).gt.ptjetmin.and.
     $        abs(eta(cshift(pjet_all(:,ijet),-1))).lt.absEtaMax) then
      !     if(ptj_all(ijet).gt.ptalljetmin.and.
      ! $        abs(yj_all(ijet)).lt.yjetmax) then
               j=j+1
               pjet(:,j)=pjet_all(:,ijet)
               ptj(j) = ptj_all(ijet)
               yj(j) = yj_all(ijet)
               phij(j) = phi_all(ijet)
          endif
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
  
      logical function is_charmed(pid)
         implicit none
         integer :: pid
         integer,dimension(32),parameter :: cpid=(/411,10411,413,10413,
     1   20413,415,431,10431,433,10433,20433,435,4122,4222,4212,4224,
     2   4214,4232,4322,4323,4412,4422,4414,4424,4432,4434,4444,5242,
     3   5422,5424,5442,5444/)
         is_charmed=.FALSE.
         if(any(cpid.eq.abs(pid)))is_charmed=.TRUE.
      end function is_charmed

      logical function charged(pid)
         implicit none
         integer :: pid
         integer,dimension(122),parameter :: cpid=(/1,2,3,4,5,6,11,13,15,
     1 24,
     1    211,9000211,100211,10211,9010211,213,10213,20213,9000213,100,213,
     1    9010213,9040213,215,10215,9000215,9010215,217,9000217,9010217,219,
     1    321,9000321,10321,100321,9010321,9020321,323,10323,20323,100323,
     1    9000323,30323,325,9000325,10325,20325,9010325,9020325,327,9010327,
     1    329,9000329,
     1    411,10411,413,10413,20413,415,431,10431,433,10433,10433,20433,435,
     1    521,10521,523,10523,20523,525,541,543,10541,10543,20543,545
     1    2212,2224,2214,1114,3222,3112,3224,3114,3312,3314,3334,4122,4222,4212,4224,
     1    4214,4232,4322,4324,4412,4422,4414,4424,4432,4434,4444,5112,5222,5114,5224,
     1    5132,5312,5314,5332,5334,5242,5422,5424,5442,5444,5512,5514,5532,5534,
     1    5554/)

         charged=.FALSE.
         if(any(cpid.eq.abs(pid)))charged=.TRUE.
      end function charged

      subroutine getdphi(p1,p2,dphi,phi1,phi2)
      implicit none
      include 'pwhg_math.h'
      real(kind=8),intent(in)   :: p1(4),p2(4)
      real(kind=8),intent(out)  :: dphi
      real(kind=8),intent(out),optional :: phi1,phi2
      real(kind=8) :: phiA, phiB,pt1,pt2

      pt1=sqrt(p1(1)**2+p1(2)**2)
      pt2=sqrt(p2(1)**2+p2(2)**2)

      if(p1(2).ge.0)then
         phiA = dacos(p1(1)/pt1)
      else
         phiA=2*pi-dacos(p1(1)/pt1)
      end if
      if(p2(2).ge.0) then
         phiB = dacos(p2(1)/pt2)
      else
         phiB = 2*pi - dacos(p2(1)/pt2)
      end if
      dphi=abs(phiA - phiB)
      if(dphi.gt.pi) dphi=2*pi-dphi
      if(present(phi1).and.present(phi2))then
         phi1=phiA
         phi2=phiB
      end if
      end subroutine getdphi
