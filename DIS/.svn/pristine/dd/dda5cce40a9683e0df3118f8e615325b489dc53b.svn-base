      program main_pythia8
      implicit none
      include "nlegborn.h"
      include "pwhg_rwl.h"
      include 'LesHouches.h'
      include "pwhg_physpar.h"
      include "pwhg_st.h"
      include 'hepevt.h'
      include "pwhg_rad.h"
      include "PhysPars.h"
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,k,l,m,iret
      integer maxev
      common/mcmaxev/maxev
      integer tune,had,hepmc,MPI,shower, QED, LO
      common/cpy8pars/tune,had,hepmc,MPI, shower, QED, LO

      real * 8 powheginput,scalupfac
      external powheginput
      integer iun
      real * 8 ub_btilde_corr,ub_remn_corr,ub_reg_corr,ub_corr
      logical NaNInEvent
      integer i1, i2

c     Pythia pp tune
      tune  = powheginput("#py8tune")

!hepmc = powheginput("#py8hepmc")
      hepmc = -1

c     Hadronization: 0 = Off, 1 = On, 2 = On with heavy hadron decays off [Default]
      had   = powheginput("#py8had")
      if(had < 0 .or. had >2) had = 2

      LO = powheginput("#LOevents")
      if (LO .ne. 1) LO = 0

C     Underlying event: 0 = Off, 1 = On [Default]
      mpi   = powheginput("#py8mpi")
      if(mpi < 0) mpi = 1          
      
c     Showers: 0 = simple dipole shower [unreccomended],
c              1 = fully local dipole shower   
c              2 = antenna shower [Vincia]
c              3 = antenna shower with "global" recoil for IF [Vincia]
      shower = powheginput("#py8shower")
      if(shower<0 .or. shower >3) shower = 1
 
c     QED shower: 0 = Off,  1= On [Default], 2= On, ``more coherent'' model used (Available only in Vicia)
      QED   = powheginput("#py8QED")

     
      if( (QED < 0) .or. (shower>2 .and. QED >2) .or.
     $     (shower <2 .and. QED > 1) ) QED = 1;
   
     
      WHCPRG='PYTHIA'
      ub_btilde_corr = powheginput('#ub_btilde_corr')
      if (ub_btilde_corr < 0d0) then
         ub_btilde_corr = 1d0
      endif
      ub_remn_corr = powheginput('#ub_remn_corr')
      if (ub_remn_corr < 0d0) then
         ub_remn_corr = 1d0
      endif
      ub_reg_corr = powheginput('#ub_reg_corr')
      if (ub_reg_corr < 0d0) then
         ub_reg_corr = 1d0
      endif
      
c particle masses for reshuffling
      st_nlight = 5
      physpar_ml(1) = 0.511d-3
      physpar_ml(2) = 0.1057d0
      physpar_ml(3) = 1.777d0
      physpar_mq(1) = 0.33d0     ! up
      physpar_mq(2) = 0.33d0     ! down
      physpar_mq(3) = 0.50d0     ! strange
      physpar_mq(4) = 1.50d0     ! charm
      physpar_mq(5) = 4.75d0   ! bottom


      call opencountunit(maxev,iun)

      call lhefreadhdr(iun)
      
      call init_hist

      call lhefinitemasses

      call pythia_init
     
      nevhep=0

      do l=1,maxev

         call lhefreadev(iun)

         NanInEvent=.false.
         do i1 =1, maxnup
            do i2=1,5
               NanInEvent=NanInEvent .or. (.not.(pup(i2,i1)==pup(i2,i1)))
            enddo
         enddo
 
         if(NanInEvent) cycle

c rescale the weight of the event depending on the rad_type (1..btilde, 2..remn, 3..reg)
c   using the ub_..._corr factors
         if (rad_type == 1) then
            ub_corr = ub_btilde_corr
         else if (rad_type == 2) then
            ub_corr = ub_remn_corr
         else if (rad_type == 3) then
            ub_corr = ub_reg_corr
         else
            ub_corr = 1d0
         endif
         if (rwl_num_weights.gt.0) then
            rwl_weights(1:rwl_num_weights)=
     $           ub_corr * rwl_weights(1:rwl_num_weights)
         endif
         xwgtup = ub_corr * xwgtup
      
           
         
     
         
        


         do m=1,5
c Insist to shower this event;
            call pythia_next(iret)
            if(iret.ne.1) then
               write(*,*) ' return code ',iret
               if(m.eq.1) then
                  write(*,*) ' Pythia could not shower this event'
                  call printleshouches
               endif
               write(*,*) ' retry'
            else
               exit
            endif
         enddo

         if(iret.eq.1) then
            call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
     1           jdahep,phep,vhep)
            if(nevhep.lt.6) then
               do j=1,nhep
                  write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
     1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
               enddo
            endif
            call pyanal
            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
               call pyaend
            endif
            if(nevhep.gt.0.and.mod(nevhep,1000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
            endif
         endif
      enddo

      write(*,*) 'At the end NEVHEP is ',nevhep

      call pythia_stat

!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
c---user's terminal calculations
      call pyaend
 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
      
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup
      end

      
      subroutine pyaend
      character * 6 vetoname
      character * 20 pwgprefix
      character * 100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto

      if(analysis_jetveto) then
         vetoname=trim(adjustl('_veto'))
      else
         vetoname=trim(adjustl(''))
      endif
      if(rnd_cwhichseed.ne.'none') then
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output-'
     1        //rnd_cwhichseed//vetoname
      else
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output'
     $        //vetoname
      endif

      call pwhgsetout
      call pwhgtopout(filename)
      end

      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf.
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      write(nlf,'(a)')'</event>'
 210  format(1p,2(1x,i8),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end

