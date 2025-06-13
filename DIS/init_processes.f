c     
c     initialization routine for ep->ej or vp->vj
c     
      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
      include 'pwhg_st.h'
      logical debug
      parameter (debug=.true.)
      integer j,i,ii,jj,k
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      integer lcharge3(-16:16)  !3*lepton charge (e-,v)
      data lcharge3 /0,3,0,3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,-3,0,-3,0/
      logical condition
      integer ferm_charge(5)
      character *3 flav(-5:5)
      data (flav(i),i=-5,5) 
     #/'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      integer max_flav
      logical emit_Wp_upper,emit_Wm_upper,emit_Wp_lower,emit_Wm_lower
      integer flst_nborn_WW,flst_nreal_WW
      integer flst_nreal_ZZ,flst_nborn_ZZ
      integer flst_nreal_WW_qq,flst_nreal_WW_qg,flst_nreal_WW_gq
      integer flst_nreal_ZZ_qq,flst_nreal_ZZ_qg,flst_nreal_ZZ_gq
      logical CKM_diag
      integer flst_born_tmp(nlegbornexternal,maxprocborn)
      logical flavequiv
      external flavequiv
      logical tag,newtag
      real * 8 powheginput
      external powheginput

      logical nc_include,cc_include
      integer channel_type  
      logical cc_only
      integer ftype,lflav
      integer id_beam
      logical, save:: allow_debug = .false.
c     
c******************************************************
c     flavor channels to be considered for 
c     e- q  -> e- q or e- q -> ve q (charged lepton beam) or 
c     v q  -> v q or v q -> e q (neutrino beam)
c      
c     sub-processes and real gluon emission (and crossings):
c     
c     NC:
c     eueu-type(channel_type=1)
c     eded-type(channel_type=2)
c     
c     all NC type contributions: (channel_type=4)
c     
c     CC:
c     euvd-type(channel_type=3)
c     
c     all CC type contributions: (channel_type=3)
!     To handle electron we fake a lepton pdf, hence the need to adjust the pdf array to go abs(ih1). 
      pdf_nparton = abs(powheginput("ih1"))
      lflav = powheginput("ih1")
      print*, 'pdf_nparton', pdf_nparton
      channel_type = 0          ! change via input file to select NC or CC

      id_beam=1             !charged lepton beam
      if ((abs(lflav).eq.12).or.(abs(lflav).eq.14).or.
     &    (abs(lflav).eq.16))   id_beam=0  !v beam
 
  
      flg_dis_isr = .true.
      flg_dis_fsr = .true.

      
      
      if(powheginput('channel_type').gt.0) then
         channel_type =powheginput('channel_type')
      endif
      


      if (channel_type.eq.0) then
         cc_include = .true.
         nc_include = .true.
         write(*,*) 
         write(*,*) ' all flavor channels are summed'   
      elseif (channel_type.eq.4) then
         cc_include = .false.
         nc_include = .true.
         write(*,*) 
         write(*,*) ' all neutral current channels are summed'   
         write(*,*) ' (no charged current channels)' 
      elseif (channel_type.eq.3) then
         cc_include = .true.
         nc_include = .false.
         write(*,*) 
         write(*,*) ' all charged current channels are considered'   
         write(*,*) ' (no neutral current channels)'   
      elseif (channel_type.ge.1.and.channel_type.le.2) then
         nc_include = .true.
         cc_include = .false.
         write(*,*) 
         write(*,*) ' one neutral-current channel considered only'  
         write(*,*) ' channel_type = ',channel_type
      else 
         write(*,*) ' '   
         write(*,*) 'ERROR: no valid entry for channel_type;'
         write(*,*) 'check your settings in powheg.input'
         stop
      endif  

      if(.not. allow_debug) then
         if(channel_type .lt. 3d0 .or. channel_type .gt. 4d0) then
            write(*,*) 'ERROR: no valid entry for channel_type;'
            write(*,*) 'only channel_type = 3 or 4 is allowed;'
            write(*,*) 'check your settings in powheg.input'
            stop
         endif  
      endif  


c*********************************************************
c     
      tag = .true.
      newtag = .true.

      if (.not.tag) then
         do i=1,nlegbornexternal
            do j=1,maxprocborn
               flst_borntags(i,j)=0
            enddo
         enddo
         do i=1,nlegbornexternal
            do j=1,maxprocreal
               flst_realtags(i,j)=0
            enddo
         enddo
      endif

      if (.not.tag) then
         write(*,*) 'jet tagging '//
     #'must be switched on'
         stop   
      endif

c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
c     BJ: attention here!       
!     flst_lightpart=4
      max_flav =  5             ! no b in initial state 
c     init:
      flst_nborn=0
      flst_nreal=0
      
      flst_nborn_WW = 0
      flst_nreal_WW = 0
      flst_nreal_ZZ = 0

      ftype = 0

      condition=.false.

c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC           BORN-TYPE GRAPHS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     consider a diagonal CKM matrix. CKM matrix elements will be added      
c     when the event is written on the Les Houches file
      CKM_diag = .true.

c charged current e q -> v q or v q -> e q:
      if (cc_include) then
      
      emit_Wp_upper = .true.
      emit_Wm_upper = .true.
      emit_Wp_lower = .true.
      emit_Wm_lower = .true.
      
c     CC case
c     e q -> v q or v q -> e q

      if (id_beam.eq.1) then !charged lepton beam
         i = lflav              !11
         ii = sign(abs(lflav) + 1,lflav) ! 12
      else ! v beam    
         i = lflav              !12
         ii = sign(abs(lflav) - 1,lflav) ! 11
      endif 
       
         do j=-max_flav,max_flav
               do jj=-max_flav,max_flav 
c              
                  if (.not.((i.eq.0) .or.(j.eq.0).or.
     #                      (ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                     ferm_charge(1) = lcharge3(i)
                     ferm_charge(2) = charge3(j)
                     ferm_charge(3) = lcharge3(ii)
                     ferm_charge(4) = charge3(jj)

                     if (CKM_diag) then
                        emit_Wp_upper = (ii.eq.i-1)
                        emit_Wm_upper = (ii.eq.i+1)
                        emit_Wp_lower = (jj.eq.j-1)
                        emit_Wm_lower = (jj.eq.j+1)                 
                     endif

!                     print*, 'ferm_charge(1)', ferm_charge(1)
!                     print*, 'ferm_charge(2)', ferm_charge(2)
!                     print*, 'ferm_charge(3)', ferm_charge(3)
!                     print*, 'ferm_charge(4)', ferm_charge(4)
!                     print*, 'emit_Wp_upper', emit_Wp_upper
!                     print*, 'emit_Wm_upper', emit_Wm_upper
!                     print*, 'emit_Wp_lower', emit_Wp_lower
!                     print*, 'emit_Wm_lower', emit_Wm_lower

                     
                     condition = 
c     W+ emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)+3)
     #                    .eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)-3)
     #                    .eq.0).and.(emit_Wm_lower)))
     #                    .or.
c     W- emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)-3)
     #                    .eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)+3)
     #                    .eq.0).and.(emit_Wp_lower)))
!                     print*, 'condition', condition
                     !stop
                     if (.not.condition) goto 600

                     if (channel_type.eq.0.or.channel_type.eq.3) then
                        condition=.true.
                     else   
                        condition = .false.                         
                     endif      !channel_type
                     
                     if (condition) then
                        flst_nborn=flst_nborn+1
                        if(flst_nborn.gt.maxprocborn) goto 999
                        flst_born(1,flst_nborn)=i
                        flst_born(2,flst_nborn)=j
                        flst_born(3,flst_nborn)=ii
                        flst_born(4,flst_nborn)=jj
                        if (tag) then
                           flst_borntags(1,flst_nborn)=1
                           flst_borntags(2,flst_nborn)=2
                           flst_borntags(3,flst_nborn)=3
                           flst_borntags(4,flst_nborn)=4
                           if (newtag) then
                              flst_borntags(1,flst_nborn)=1
                              flst_borntags(2,flst_nborn)=2
                              flst_borntags(3,flst_nborn)=1 
                              flst_borntags(4,flst_nborn)=2 
                           endif
                        endif
                     endif
                  endif
 600              continue
c               enddo
            enddo
c     enddo
         enddo
         flst_nborn_WW = flst_nborn

         if (debug) then
            write(*,*) ' Born processes: CC ',flst_nborn
            do j=1,flst_nborn
               write(*,*) 'proc #',j,':',(flst_born(k,j),k=1,nlegbornexternal)
            enddo
         endif

      endif                     !cc_include
      if (.not.nc_include) goto 211


c     neutral current:
      if (nc_include) then      !e q -> e q 

         i = lflav              !11        
c     do i=-max_flav,max_flav
         do j=-max_flav,max_flav
c     ~          do j=1,1
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = lcharge3(i)
               ferm_charge(2) = charge3(j)
               ferm_charge(3) = lcharge3(i)
               ferm_charge(4) = charge3(j)

               if (channel_type.eq.0.or.channel_type.eq.4) then
                  condition=.true.
               else     
c     use value of i/j to select appropriate flavor channel    
                  ftype = 2-mod(abs(j),2)         
                  k = -ftype+3

                  if (k.eq.channel_type) then
                     condition=.true.
                  else
                     condition = .false.
                  endif         !k
                  
               endif            !channel_type

               if (condition) then                      
                  flst_nborn=flst_nborn+1
                  
                  if(flst_nborn.gt.maxprocborn) goto 999
                  flst_born(1,flst_nborn)=i
                  flst_born(2,flst_nborn)=j
                  flst_born(3,flst_nborn)=i
                  flst_born(4,flst_nborn)=j
                  if (tag) then
                     flst_borntags(1,flst_nborn)=1
                     flst_borntags(2,flst_nborn)=2
                     flst_borntags(3,flst_nborn)=3
                     flst_borntags(4,flst_nborn)=4
                     if (newtag) then
                        flst_borntags(1,flst_nborn)=1
                        flst_borntags(2,flst_nborn)=2
                        flst_borntags(3,flst_nborn)=1 
                        flst_borntags(4,flst_nborn)=2 
                     endif      !newtag
                  endif         !tag
               endif            !condition

            endif               ! no gluons
c     enddo
         enddo
         flst_nborn_ZZ = flst_nborn             

         if (debug) then
            write(*,*) ' Born processes: NC ',flst_nborn-flst_nborn_WW
            do j=flst_nborn_WW+1,flst_nborn
               write(*,*) 'proc #',j-flst_nborn_WW,':',
     #(flst_born(k,j),k=1,nlegbornexternal)
            enddo
         endif

      endif                     !nc_include

! Number of particles, including resonances
      flst_bornlength = nlegbornexternal
! Number of final-state external particles 
      flst_numfinal   = nlegbornexternal-2
! resonance assigmnents, 0=all come from production
      flst_bornres    = 0

 211  continue

c     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC  REAL EMISSION GRAPHS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     consider a diagonal CKM matrix. CKM matrix elements will be added      
c     when the event is written on the Les Houches file
      CKM_diag = .true.

c     charged current:
      if (cc_include) then
         
         emit_Wp_upper = .true.
         emit_Wm_upper = .true.
         emit_Wp_lower = .true.
         emit_Wm_lower = .true.


c     charged current of type e q -> v q g or v q -> e q g:  
      if (id_beam.eq.1) then !charged lepton beam
         i = lflav              !11
         ii = sign(abs(lflav) + 1,lflav) ! 12
      else ! v beam    
         i = lflav              !12
         ii = sign(abs(lflav) - 1,lflav) ! 11
      endif 

         do j=-max_flav,max_flav
            do jj=-max_flav,max_flav               
c     
               if (.not.((i.eq.0) .or.(j.eq.0).or.
     #(ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon

                  ferm_charge(1) = lcharge3(i)
                  ferm_charge(2) = charge3(j)
                  ferm_charge(3) = lcharge3(ii)
                  ferm_charge(4) = charge3(jj)

                  if (CKM_diag) then
                     emit_Wp_upper = (ii.eq.i-1)
                     emit_Wm_upper = (ii.eq.i+1)
                     emit_Wp_lower = (jj.eq.j-1)
                     emit_Wm_lower = (jj.eq.j+1)                 
                  endif
                  
                  condition = 
c     W+ emission from upper leg                        
     #(((ferm_charge(1)-(ferm_charge(3)+3)
     #.eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #((ferm_charge(2)-(ferm_charge(4)-3)
     #.eq.0).and.(emit_Wm_lower)))
     #.or.
c     W- emission from upper leg                        
     #(((ferm_charge(1)-(ferm_charge(3)-3)
     #.eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #((ferm_charge(2)-(ferm_charge(4)+3)
     #.eq.0).and.(emit_Wp_lower)))

                  if (.not.condition) goto 800

                  if (channel_type.eq.0.or.channel_type.eq.3) then
                     condition=.true.
                  else      
                     condition = .false.
                  endif         !channel_type
                  
                  if (condition) then
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     flst_real(1,flst_nreal)=i
                     flst_real(2,flst_nreal)=j
                     flst_real(3,flst_nreal)=ii
                     flst_real(4,flst_nreal)=jj
                     flst_real(5,flst_nreal)=0 ! gluon
                     if (tag) then
                        flst_realtags(1,flst_nreal)=1
                        flst_realtags(2,flst_nreal)=2
                        flst_realtags(3,flst_nreal)=3
                        flst_realtags(4,flst_nreal)=4
                        flst_realtags(5,flst_nreal)=0
                        if (newtag) then
                           flst_realtags(1,flst_nreal)=1
                           flst_realtags(2,flst_nreal)=2
                           flst_realtags(3,flst_nreal)=1 
                           flst_realtags(4,flst_nreal)=2 
                           flst_realtags(5,flst_nreal)=2
                        endif
                     endif
                  endif
               endif
 800           continue
c     enddo
            enddo
c     enddo
         enddo
         flst_nreal_WW_qq = flst_nreal
         
         
ccccccc

c     e g -> v q q or v g -> e q q:
      if (id_beam.eq.1) then !charged lepton beam
         i = lflav              !11
         ii = sign(abs(lflav) + 1,lflav) ! 12
      else ! v beam    
         i = lflav              !12
         ii = sign(abs(lflav) - 1,lflav) ! 11
      endif        
         
c     loop on only HALF of the incoming lower-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.
         do j=-1, -max_flav, -1             
            do jj = 1 , max_flav
               if (.not.((i.eq.0).or.(j.eq.0).or.
     #(ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                  ferm_charge(1) = lcharge3(i)
                  ferm_charge(2) = 0
                  ferm_charge(3) = lcharge3(ii)
                  ferm_charge(4) = charge3(jj)
                  ferm_charge(5) = charge3(j)
                  
                  condition = ((ferm_charge(3)+ferm_charge(4)+ferm_charge(5)
     $                 -ferm_charge(1)-ferm_charge(2)).eq.0d0)
                  
  
                  
                  if (.not.condition) goto 803

                  if(CKM_diag) then
                     if(abs(j) == 1 .and. abs(jj) .ne. 2) goto 803
                     if(abs(j) == 2 .and. abs(jj) .ne. 1) goto 803
                     if(abs(j) == 3 .and. abs(jj) .ne. 4) goto 803
                     if(abs(j) == 4 .and. abs(jj) .ne. 3) goto 803
                     if(abs(j) == 5 .and. abs(jj) .ne. 6) goto 803
                     if(abs(j) == 6 .and. abs(jj) .ne. 5) goto 803
                  endif

c     
                  if (channel_type.eq.0.or.channel_type.eq.3) then
                     condition=.true.
                  else     
                     condition = .false.
                  endif         !channel_type

                  if (condition) then
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     flst_real(1,flst_nreal)=i
                     flst_real(2,flst_nreal)=0 ! gluon
                     flst_real(3,flst_nreal)=ii
                     flst_real(4,flst_nreal)=jj
                     flst_real(5,flst_nreal)=j
                     if (tag) then
                        flst_realtags(1,flst_nreal)=1
                        flst_realtags(2,flst_nreal)=0
                        flst_realtags(3,flst_nreal)=3
                        flst_realtags(4,flst_nreal)=4
                        flst_realtags(5,flst_nreal)=2
                        if (newtag) then
                           flst_realtags(1,flst_nreal)=1
                           flst_realtags(2,flst_nreal)=2
                           flst_realtags(3,flst_nreal)=1 
                           flst_realtags(4,flst_nreal)=2 
                           flst_realtags(5,flst_nreal)=2
                        endif                           
                     endif

                  endif
               endif
 803           continue   
            enddo
c     enddo
         enddo
c     enddo

         flst_nreal_WW_gq = flst_nreal
         flst_nreal_WW = flst_nreal
         
         if (debug) then
            write(*,*) ' real processes: CC ',flst_nreal
            do j=1,flst_nreal
               write(*,*) 'proc #',j,':',(flst_real(k,j),k=1,nlegreal)
            enddo
         endif

      endif                     !cc_include
      if (.not.nc_include) goto 111

cccccccc
c     
c     neutral current:
      if (nc_include) then

c     e q -> e q g
         i = lflav              !11   
c     do i=-max_flav,max_flav
         do j=-max_flav,max_flav
c     ~          do j=1,1

            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = lcharge3(i)
               ferm_charge(2) = charge3(j)
               ferm_charge(3) = lcharge3(i)
               ferm_charge(4) = charge3(j)

               if (channel_type.eq.0.or.channel_type.eq.4) then
                  condition=.true.
               else     
c     use value of i/j to select appropriate flavor channel    
                  ftype = 2-mod(abs(j),2)         
                  k = -ftype+3

                  if (k.eq.channel_type) then
                     condition=.true.
                  else
                     condition = .false.
                  endif         !k
                  
               endif            !channel_type

               if (condition) then                      
                  flst_nreal=flst_nreal+1
                  
                  if(flst_nreal.gt.maxprocreal) goto 998
                  flst_real(1,flst_nreal)=i
                  flst_real(2,flst_nreal)=j
                  flst_real(3,flst_nreal)=i
                  flst_real(4,flst_nreal)=j
                  flst_real(5,flst_nreal)=0 ! gluon
                  if (tag) then
                     flst_realtags(1,flst_nreal)=1
                     flst_realtags(2,flst_nreal)=2
                     flst_realtags(3,flst_nreal)=3
                     flst_realtags(4,flst_nreal)=4
                     flst_realtags(5,flst_nreal)=0
                     if (newtag) then
                        flst_realtags(1,flst_nreal)=1
                        flst_realtags(2,flst_nreal)=2
                        flst_realtags(3,flst_nreal)=1 
                        flst_realtags(4,flst_nreal)=2 
                        flst_realtags(5,flst_nreal)=2
                     endif      !newtag
                  endif         !tag
               endif            !condition

            endif               ! no gluons
c     enddo
         enddo
         flst_nreal_ZZ_qq = flst_nreal
cccccc
c     e g -> e q q
c     loop on only HALF of the incoming lower-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.

         i = lflav              !11
c     do i=-max_flav,max_flav
         do j=1,max_flav  
c     ~          do j=1,-1  
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = lcharge3(i)
               ferm_charge(2) = 0
               ferm_charge(3) = lcharge3(i)
               ferm_charge(4) = charge3(j)
               ferm_charge(5) = -charge3(j)

               if (channel_type.eq.0.or.channel_type.eq.4) then
                  condition=.true.
               else     
                  
c     use value of i/j to select appropriate flavor channel    
                  ftype = 2-mod(abs(j),2)              
                  k = -ftype+3

                  if (k.eq.channel_type) then
                     condition=.true. 
                  else  
                     condition = .false.
                  endif         !k   
               endif            !channel_type


               if (condition) then                
                  flst_nreal=flst_nreal+1

                  if(flst_nreal.gt.maxprocreal) goto 998
                  flst_real(1,flst_nreal)=i
                  flst_real(2,flst_nreal)=0 ! gluon
                  flst_real(3,flst_nreal)=i
                  flst_real(4,flst_nreal)=j
                  flst_real(5,flst_nreal)=-j
                  if (tag) then
                     flst_realtags(1,flst_nreal)=1
                     flst_realtags(2,flst_nreal)=0
                     flst_realtags(3,flst_nreal)=3
                     flst_realtags(4,flst_nreal)=4
                     flst_realtags(5,flst_nreal)=2
                     if (newtag) then
                        flst_realtags(1,flst_nreal)=1
                        flst_realtags(2,flst_nreal)=2
                        flst_realtags(3,flst_nreal)=1 
                        flst_realtags(4,flst_nreal)=2 
                        flst_realtags(5,flst_nreal)=2
                     endif
                  endif
               endif
            endif
         enddo
c     enddo 


         flst_nreal_ZZ_gq = flst_nreal

         if (debug) then
            write(*,*) ' real processes: NC ',flst_nreal-flst_nreal_WW
            do j=flst_nreal_WW+1,flst_nreal
               write(*,*) 'proc #',j-flst_nreal_WW,':',
     #(flst_real(k,j),k=1,nlegreal)
            enddo
         endif

      endif                     !nc_include

!Number of external particles, including resonces
      flst_reallength = nlegbornexternal+1
!all comes from production
      flst_realres    = 0

 111  continue
c     
      
      
      call buildresgroups(flst_nborn,nlegborn,flst_bornlength,
     1     flst_born,flst_bornres,flst_bornresgroup,flst_nbornresgroup)
      call buildresgroups(flst_nreal,nlegreal,flst_reallength,
     1     flst_real,flst_realres,flst_realresgroup,flst_nrealresgroup)

      st_nlight=5
!   call init_couplings    [already in init_phys.f]

      return
 998  write(*,*) 'init_processes: increase maxprocreal:',maxprocreal
      stop
 999  write(*,*) 'init_processes: increase maxprocborn:',maxprocborn
      stop
      
      end



