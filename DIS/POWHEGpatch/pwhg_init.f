      subroutine pwhginit
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      real * 8 powheginput
      character * 3 whichpdfpk
      external powheginput,whichpdfpk
      integer i1,n1,n2
      call init_flstarrays
c whether the remnant use the same grid as btilde
      flg_samexgridasbtilde = .false.
c whether the kinematics for massive emitters is partly done
c in real * 16
      flg_mvemquadprec = .false.
c     flg_newrad is no longer used. It should be deleted from the include file too
c      flg_newrad = .not.powheginput("#newrad").eq.0
      flg_allrad = powheginput("#allrad").eq.1
      flg_debug=.false.
c default for scheme in Virtual
c flg_drscheme=.false. Traditional MSbar
c flg_drscheme=.true.  Dimensional reduction
      flg_drscheme=.false.
c Only photon radiation. NO QCD
      flg_QEDonly=.false.
c On if there are electromagnetic corrections
      flg_with_em=.false.
c default is 6; if em is included, it is 22 to include photons
      pdf_nparton=6
      if(powheginput("#flg_debug").eq.1) flg_debug=.true.
c whether to output negative weights or not
      flg_withnegweights=.true.
      if(powheginput("#withnegweights").eq.0)flg_withnegweights=.false.
c See if we have weighted events
      flg_weightedev=.false.
      if(powheginput("#bornsuppfact").gt.0) flg_weightedev=.true.
      if(powheginput("#ptsupp").gt.0) then
         write(*,*) ' ******** WARNING: ptsupp is deprecated'
         write(*,*) ' ******** Replace it with bornsuppfact'
         call flush(6)
         call exit(-1)
      endif
c if true, suppress emitted particle energy with respect to the emitter
c in all cases (i.e. not only for g g) in FSR. Thus, consider also
c the emitter-emitted pairs qbar q, g q. 
      flg_doublefsr=.false.
      if(powheginput("#doublefsr").gt.0) flg_doublefsr=.true.
c     Set to true to remember and use identical values of the computed 
c     amplitudes, for Born, real and virtual contributions
      flg_smartsig=.true.
      flg_bornsmartsig=.true.
c     By default the cache machinary discards amplitudes which are zero
c     for a sufficient number of calls (see cache_similar.f, around line 180)
c     set to .false. in your init_processes if there are hard cuts in the
c     phase space (setting jakobian to zero) 
      flg_cachediscardzeroamps=.true.
c the following must be set to true while gen is called
      flg_ingen = .false.
      if(powheginput("#smartsig").eq.0) then
         flg_smartsig=.false.
         flg_bornsmartsig=.false.
      endif
      if(powheginput("#bornsmartsig").eq.0) flg_bornsmartsig=.false.
c     if true also counterterm are output in NLO tests (default)
      flg_withsubtr=.true.
      if(powheginput("#withsubtr").eq.0) flg_withsubtr=.false.
c     The following turns on mechanism to deal with Born zeroes
      flg_bornzerodamp=.false.
c If the damp function has to be called, the following must be on
      flg_withdamp=.false.
c Traditional mechanism: either withdamp 1 or hfact >0 turn on
c remnant calculation. This mechanism is prone to error, since one may
c want hfact without bornzerodamp. We keep it for backward compatibility.
c In the future, we should use the flags listed below
      if(powheginput("#hfact").gt.0d0) then
         flg_withdamp=.true.
         flg_bornzerodamp=.true.
      endif
      if(powheginput("#withdamp").eq.1) then
         flg_withdamp=.true.
         flg_bornzerodamp=.true.
      endif
c Modern mechanism;
      if(powheginput("#bornzerodamp").eq.1) then
         flg_withdamp=.true.
         flg_bornzerodamp=.true.
      endif
c bornzerodamp can be forced to 0 anyhow, if required explicitly
      if(powheginput("#bornzerodamp").eq.0) then
         flg_bornzerodamp=.false.
      endif
c hdamp will replace the old hfact in new implementations
      if(powheginput("#hdamp").gt.0) then
         flg_withdamp=.true.
      endif
c     If set do only the Born term
      flg_bornonly=.false.
      if (powheginput("#bornonly").eq.1) flg_bornonly=.true.
      flg_virtonly=.false.
      if (powheginput("#virtonly").eq.1) flg_virtonly=.true.
      flg_softonly=.false.
      if (powheginput("#softonly").eq.1) flg_softonly=.true.
      flg_LOevents=.false.
      if (powheginput("#LOevents").eq.1) then
         flg_LOevents=.true.
         flg_bornonly=.true.
      endif
      flg_novirtual=.false.
      if (powheginput("#novirtual").eq.1) then
         if(flg_virtonly) then
            write(*,*) ' pwhginit: you have specified virtonly'
            write(*,*) '           and novirtual at the same time'
            write(*,*) ' cannot proceed, exiting ...'
            call exit(-1)
         endif
         if(flg_softonly) then
            write(*,*) ' pwhginit: you have specified virtonly'
            write(*,*) '           and softonly at the same time'
            write(*,*) ' cannot proceed, exiting ...'
            call exit(-1)
         endif
         flg_novirtual=.true.
      endif
c     for_reweighting gives settings that are supposed to be best if we first
c     generate events without the virtual and reweight subsequently.
      if (powheginput("#for_reweighting") == 1) then
         if(powheginput('#compute_rwgt') == 1) then
            write(*,*) ' you have specified compute_rwgt'//
     1           ' and for_reweighting at the same time;'
            write(*,*) ' We ignore  for_reweighting.'
            flg_for_reweighting = .false.
         else
            if(powheginput("#virtonly") == -1d6) then
               write(*,*) 'You have specified for_reweighting:'
               write(*,*) 'we also set flg_novirtual = .true.'
            endif
            flg_novirtual = .true.
            flg_for_reweighting = .true.
         endif
      else
         flg_for_reweighting = .false.
      endif
c     initialize random number sequence
      i1=powheginput('#iseed')
      if (i1.lt.0) i1=0
      n1=powheginput('#rand1')
      if (n1.lt.0) n1=0
      n2=powheginput('#rand2')
      if (n2.lt.0) n2=0
c select which upper bounding function form
      rad_iupperisr=powheginput("#iupperisr")
      if(rad_iupperisr.lt.0) rad_iupperisr=1
      rad_iupperfsr=powheginput("#iupperfsr")
      if(rad_iupperfsr.lt.0) rad_iupperfsr=2
c info on pdf for each event
      flg_pdfreweight=.false.
      if(powheginput("#pdfreweight").eq.1)flg_pdfreweight=.true.
c infos to write LH file with infos to reweight events
      flg_reweight=.not.(powheginput('#storeinfo_rwgt').eq.0)
c infos to read LH file with infos to reweight events,
c and produce a new LH file with updated weights
      flg_newweight=powheginput('#compute_rwgt').eq.1
c Is MiNLO required?
      flg_minlo=.false.
      st_bornorder=0
      if(powheginput("#minlo").eq.1) then
         flg_minlo=.true.
c set st_bornorder to an impossible value; later on we check
c if the user has set it. This is only used by minlo
         st_bornorder=-1000
      endif

      pdf_alphas_from_PDF = powheginput("#alphas_from_pdf") == 1
      if(whichpdfpk() == 'mlm') then
         if(pdf_alphas_from_PDF) then
            write(*,*) 'pwhg_init: ERROR!'
            write(*,*) 'Error: we are using mlm PDF'
            write(*,*) 'and and alphas_from_PDF is set to 1,'
            write(*,*) 'not allowed, exiting...'
            call exit(-1)
         endif
      endif
      if(powheginput("#alphas_from_lhapdf") == 1) then
         if(powheginput("#alphas_from_pdf") == 0) then
            write(*,*) 'pwhg_init: ERROR!'
            write(*,*) 'alphas_from_pdf=0 and alphas_from_lhapdf=1,'
            write(*,*) 'confusing powheg.input, exiting ...'
            call exit(-1)
         endif
c     Here for upper compatibility
         if(.not. whichpdfpk() == 'lha') then
            write(*,*)' error: alphas_from_lhapdf set'//
     1           ' to 1, but we are not using lhapdf;'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         write(*,*) '************************'//
     1        'WARNING **********************'
         write(*,*) 'alphas_from_lhapdf is '//
     1        'deprecated; use instead alphas_from_pdf'
         pdf_alphas_from_PDF = .true.
      endif

      
c     whether to correct for upper bound violations in the generation
c     of btilde and remnant events 
      flg_ubexcess_correct = powheginput("#ubexcess_correct") ==1           
c
c Light flavours, i.e. flavours that can undergo collinear splitting
      flst_collparticles = 1000000
      flst_ncollparticles = 11
      flst_collparticles(1:flst_ncollparticles)=
     1     (/ 0,1,2,3,4,5,11,13,15,21,22 /)
      call setrandom(i1,n1,n2)
c     assign a default id for the process at hand
c     if the user want to assign different id's
c     to each subprocess, he/she can reassign lprup(1)
c     inside the user-defined subroutine init_processes
      lprup(1)=10001
c     initialize physical parameters
      call init_phys
      if(flg_minlo.and.st_bornorder.eq.-1000) then
         write(*,*) '************************************'
         write(*,*) ' ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*) ' st_bornorder should be set by the user'
         write(*,*) ' it should equal the powers of alpha_S'
         write(*,*) ' present at the Born level'
         write(*,*) ' Exiting ...'
         call pwhg_exit(-1)
      endif
c ID of beam particles 1 and 2 (pdg convention)
c proton:
      ebmup(1)=kn_beams(0,1)
      ebmup(2)=kn_beams(0,2)
      if(abs(pdf_ih1).eq.1) then
c proton and antiproton
        idbmup(1) = 2212*pdf_ih1
      elseif(abs(pdf_ih1).eq.2 .or. abs(pdf_ih1).eq.22) then
c neutron and antineutron (2 or 22 for when input pdfs are proton or 
c neutron. Option 2 assumes proton and performs and isospin 
c transformation within lhapdf6if.cc
         idbmup(1) = 2112*pdf_ih1/abs(pdf_ih1)
c leptons
      elseif(abs(pdf_ih1).ge.11 .or. abs(pdf_ih1).le.16) then                  
         idbmup(1) = pdf_ih1
      else
         write(*,*) ' unimplemented hadron 1 ',pdf_ih1
         call pwhg_exit(-1)
      endif
      if(abs(pdf_ih2).eq.1) then
        idbmup(2) = 2212*pdf_ih2
      elseif(abs(pdf_ih2).eq.2 .or. abs(pdf_ih2).eq.22) then
         idbmup(2) = 2112*pdf_ih2/abs(pdf_ih2)
      elseif(abs(pdf_ih1).eq.11) then
c     electrons
         idbmup(1) = 11*pdf_ih1/abs(pdf_ih1)
      else
         write(*,*) ' unimplemented hadron 2 ',pdf_ih2
         call pwhg_exit(-1)
      endif
c pdf group; negative to use internal herwig pdf's for showering
      pdfgup(1)=-1
      pdfgup(2)=-1
c pdf set
      pdfsup(1)=-1
      pdfsup(2)=-1
c If either weighted events or event with negative weights are
c required, use idwtup=-4. Thuse, in both these cases, the average of the
c event weight is the total cross section.
c Otherwise the weight is set to 1 for each event.
c User processes may override these choices. In particular, using -4 in
c all cases is recommended for new processes. Using 3 is left here for
c compatibility with older implementations.
      if(flg_withnegweights.or.flg_weightedev) then
         idwtup = -4
      else
         idwtup = 3
      endif
c number of user subprocesses
c Irrelevant if idwtup=+-3,+-4
      nprup = 1
      call bbinit
c now the cross section is available
      if(flg_weightedev) then
         xsecup(1)=-1
         xerrup(1)=-1
      else
c     Provided by Mauro Chiesa, 21/05/2021
         call comp_radtot(rad_tot,rad_etot)         
         xsecup(1)=rad_tot  *rad_branching
         xerrup(1)=rad_etot *rad_branching
      endif
      xmaxup(1)=1
      end

      subroutine init_flstarrays
      include 'nlegborn.h'
      include 'pwhg_flst.h'

      flst_bornlength = nlegborn
      flst_borntags = 0
      flst_bornres = 0

      flst_reallength = nlegreal
      flst_realtags = 0
      flst_realres = 0

c not used by default
      flst_nbornresgroup = 1
      flst_bornresgroup = 1

      flst_nrealresgroup = 1
      flst_realresgroup = 1

      end


c     Provided by Mauro Chiesa, 21/05/2021
      subroutine comp_radtot(radtot,radetot)
      implicit none
      double precision radtot,radetot
      include 'pwhg_mintw.h'
      integer j
      radtot =0d0
      radetot=0d0
      do j=1,mintw_nptrs
         radtot  = radtot  + mintw_ptrs(j)%tot
         radetot = radetot + mintw_ptrs(j)%etot**2
      enddo
      radetot=sqrt(radetot)
      end subroutine comp_radtot
      
