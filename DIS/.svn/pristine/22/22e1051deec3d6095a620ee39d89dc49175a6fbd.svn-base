      subroutine testsuda
c tests the implementation of the generation of radiation
c It computes the integral of the R/B (real over born)
c times theta(pt2(r)-pt2) as a function of pt2,
c and then computes the negative exponent of the integral
c (which is the probability of no emission).
c It then computes the probability of radiation above a given pt cut.
c This is related to the no-emission probability (1-m the above).
c the two results are compared.
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_par.h'
      include 'pwhg_flg.h'
      include 'pwhg_em.h'
      integer ncalls1,ncalls2,nbins
      parameter (ncalls1=200000,ncalls2=200000,nbins=40)
      real * 8 hist1(nbins),hist2(nbins),histsq1(nbins),av1,avsq1,err1       
      real * 8 sig,born,ptsq,www,xjac,tmp
      real * 8 random,pwhg_pt2,powheginput
      real * 8 ptsqmin
      character * 50 string
      external random,pwhg_pt2,powheginput
      logical refuse_pdf
      integer j,k,firstreg,lastreg,fl1,fl2,alr
      real * 8 csi,y,azi, branchfactor
      integer iun
      logical ini
      data ini/.true./
      save ini,iun
      logical is_charged,is_coloured
      external is_charged,is_coloured
      integer em,flemitter
      integer radregsave
      save radregsave
      logical, save:: ISR = .true.
c      real * 8 wwwref

      if(ini) then         
c     if allrad is set, gen_radiation stores the kn_csi's, kn_y's, kn_azi's and 
c     only the last one would be available globally.
c     To avoid this, allrad should be false
         if (flg_allrad) then
            write(*,*) 'The testsuda check cannot run with allrad set to 1'
            write(*,*) 'Set it to 0, or comment it in the powheg.input file'
            call pwhg_exit(-1)
         endif
         radregsave=0
         call newunit(iun)
         open(unit=iun,file='testsuda.top',
     #        status='unknown')
         firstreg=powheginput("#radregion")
         if(firstreg.le.0) then
            firstreg=1
            lastreg=rad_nkinreg
         else
            lastreg=firstreg
            radregsave=firstreg
         endif
         ini=.false.
      endif
      do k=1,nbins
         hist1(k)=0
         hist2(k)=0
         histsq1(k)=0
      enddo

c$$$      if (flg_QEDonly) then
c$$$c     check photon Sudakov ONLY
c$$$         em_rad_on = .true.
c$$$      else
c$$$c     check QCD Sudakov ONLY
c$$$         em_rad_on = .false.
c$$$      endif
c$$$      em_rad_on = .true.
      
c$$$      if(.not.is_coloured(flemitter).and.is_charged(flemitter))
c$$$     1     then
c$$$         em_rad_on = .true.
c$$$      else
c$$$         em_rad_on = .false.
c$$$      endif


      if (em_rad_on) then
         ptsqmin = rad_ptsqmin_em
      else
         ptsqmin = rad_ptsqmin
      endif


      write(*,*) ' computing sudakov integral'
      call randomsave
      do rad_kinreg=firstreg,lastreg
         if(rad_kinreg_on(rad_kinreg)) then
            if(rad_kinreg.eq.1.and.ISR) then               
               kn_emitter = 0
               ! For dis kn_emitter can only be 2
               if (flg_dis_isr) kn_emitter = 2
               fl1=flst_born(1,rad_ubornsubp)
               fl2=flst_born(2,rad_ubornsubp)
               if((.not.is_coloured(fl1).and..not.is_coloured(fl2))
     1          .and.(is_charged(fl1).or.is_charged(fl2))) then
c     For processes with partons (colored)  in the initial state, we never enter here. 
c     Only when a lepton is incoming, the following apply
                  em_rad_on = .true.
               else
                  if (flg_QEDonly) then
c     in this case, ONLY real diagrams with photon radiation should be active in the 
c     init_processes.f file.     NO QCD radiation at all
                     em_rad_on = .true.
                  else
                     em_rad_on = .false.
                  endif
               endif
            else
c     final state
               em=flst_lightpart+rad_kinreg-2
               kn_emitter = em
               flemitter=flst_born(em,rad_ubornsubp)
            endif
            do j=1,ncalls1
               tmp=random()
               kn_csitilde=tmp**2
               kn_csitilde=kn_csitilde*(1-par_isrtinycsi)+par_isrtinycsi
               ! do we need to do something for csi for ISR?
               
               xjac=2*tmp
               tmp=random()
               kn_y=1-2*tmp
               xjac=xjac*2
               xjac=xjac*1.5d0*(1-kn_y**2)*(1-par_fsrtinyy)

               kn_y=1.5d0*(kn_y-kn_y**3/3)
               kn_azi=2*pi*random()
               xjac=xjac*2*pi
               kn_jacreal = 0d0
               if(rad_kinreg.eq.1.and.ISR) then
c     initial state radiation
                  call gen_real_phsp_isr_rad0
               else
c     final state radiation
                  call gen_real_phsp_fsr_rad0
               endif
               ptsq=pwhg_pt2()
               if(ptsq.gt.ptsqmin) then
                  call set_rad_scales(ptsq)
c     for photon radiation, we need to go to smaller value of pt with respect
c     to gluon radiation. So we apply refuse_pdf only to QCD radiation
c     and not to photon radiation
                  if(flg_QEDonly .or. .not.refuse_pdf()) then
                     call sigborn_rad(born)
                     call sigreal_rad(sig)
                     if(kn_jacreal.ne.0d0) then
                        www=sig/born*kn_jacreal*kn_csimax*xjac
                     else
                        www = 0d0
                     endif
!                     print*, www
!                     stop
c                     write(54,*) www
c                     read(54,*) wwwref
c                     if (abs(www/wwwref-1).gt.1d-6) then
c                        write(*,*) 'Error in testsuda, www ratio: ',www/wwwref
c                     endif
                     do k=1,nbins
                        if(log(ptsq/ptsqmin).gt.
     1                       k*log(kn_sbeams/ptsqmin)/nbins)
     2                       then
                           hist1(k)=hist1(k)+www
                           histsq1(k)=histsq1(k)+www**2
                        endif
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
      write(*,*) ' generating radiation'
      do j=1,ncalls2
         if (mod(j,1000).eq.0) then
            write(*,*) "generating radiation ",j,"/",ncalls2
         endif
         call gen_radiation
c     We should get here the generated ptsq, since the highest bid is not
c     necessarily the last one generated
         call handle_radiations('next',alr,ptsq,csi,y,azi,branchfactor)
         if (radregsave.gt.0) then
            rad_kinreg=radregsave
         endif

         if(ptsq.ge.ptsqmin) then
            do k=1,nbins
               if(log(ptsq/ptsqmin).gt.
     $              k*log(kn_sbeams/ptsqmin)/nbins) then
                  hist2(k)=hist2(k)+1
               endif
            enddo
         endif
      enddo
      call randomrestore 
      write(iun,*) ' newplot'
      write(iun,*) ' Title top "',rad_ubornsubp,'"'
      call  from_number_to_madgraph
     #         (flst_bornlength(rad_ubornsubp),flst_born(1,rad_ubornsubp),-1,string)
      write(iun,*) '(',string
      write(iun,*) ' set order x y'
      write(iun,*) ' set scale x log'
      do k=1,nbins
         av1=hist1(k)/ncalls1
         write(iun,*) exp(k*log(kn_sbeams/ptsqmin)/nbins)*
     $        ptsqmin,     exp(-av1)
      enddo
      write(iun,*) ' join'
      do k=1,nbins
         av1=hist1(k)/ncalls1
         avsq1=histsq1(k)/ncalls1
         err1=sqrt((avsq1-av1**2)/ncalls1)        
         write(iun,*) exp(k*log(kn_sbeams/ptsqmin)/nbins)*
     $        ptsqmin,     exp(-av1+err1)
      enddo
      write(iun,*) ' join'
      do k=1,nbins
         av1=hist1(k)/ncalls1
         avsq1=histsq1(k)/ncalls1
         err1=sqrt((avsq1-av1**2)/ncalls1)
         write(iun,*) exp(k*log(kn_sbeams/ptsqmin)/nbins)*
     $        ptsqmin,     exp(-av1-err1)
      enddo
      write(iun,*) ' join' 
      write(iun,*) ' set order x y dy'
      do k=1,nbins         
         write(iun,*) exp(k*log(kn_sbeams/ptsqmin)/nbins)*
     1        ptsqmin,     1-hist2(k)/ncalls2,
     2        sqrt(dble(hist2(k)))/ncalls2
      enddo
      write(iun,*) ' join dashes'
      write(iun,*) ' plot'
      call flush(iun)
!      stop
      end


