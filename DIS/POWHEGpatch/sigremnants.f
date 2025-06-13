c Function to be integrated with mint, used to generate
c non-singular contributions using standard mint-gen method
c These contributions arise from real graphs without singular regions,
c or, when damping of R/B value is adopted, from the remnant of the
c damping
      function sigremnant(iresgroup,xx,ww,ifirst,imode,retval,retval0)
c retval is the function return value
c retvavl0 is an 'avatar' function the has similar value, but is much
c easier to compute (i.e. the Born term in this case)
c imode = 0 compute retval0 only.
c imode = 1 compute retval, retval0
c return value: output, 0: success; 1: retval0 was not computed
c                 (this function does not support an avatar function)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      integer sigremnant,iresgroup,imode
      real * 8 retval,retval0,xx(ndiminteg),ww
      integer ifirst
      real * 8 xrad(3)
      real * 8 xborn(ndiminteg-3)
      integer j,alr
      real * 8 signed_sig,abs_sig,tmp
      real * 8 jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     1         jac_over_csi_s,jac_over_csi_coll
      real * 8, save ::  rrr_damp_rem_arr(maxalr),
     1     rrr_rem_arr(maxalr),rrr_damp_rem_tot
      real * 8 xjac
      logical valid_emitter
      external valid_emitter
      logical pwhg_isfinite
      external pwhg_isfinite
      flst_ibornresgroup = iresgroup
      call setup_resgroupstuff
      sigremnant = 1
      if(ifirst.eq.2) then
c     rrr_famp_rem_tot is the sum of absolute values of contributions to the
c     cross section, each with its suppression factor included.
         retval=rrr_damp_rem_tot
         if(.not.flg_ingen) then
c     This call stores, for each contribution (remnant in this case) the full array
c     of contributions, all including their suppression factor,
c     that is accumulated into sum of the contributions, sum
c     of the absolute values, etc.. The .not. flg_ingen clause checks that we are not
c     generating an event
            call adduptotals(rrr_damp_rem_arr,flst_nalr)
         else
c If we are generating an event, this call sets up arrays to            
c to pick a flavour configuration for the event.
            call storeradarray(rrr_damp_rem_arr)
         endif
         if(flg_nlotest) call pwhgaccumup
         return
      endif
      do j=1,ndiminteg-3
         xborn(j)=xx(j)
      enddo
      do j=1,3
         xrad(j)=xx(ndiminteg-3 + j)
         rad_xradremn(j)=xrad(j)
      enddo
      kn_emitter=0
      call gen_born_phsp(xborn)
c set scales
      call setscalesbtilde
c the following is needed to compute soft and collinear limits
      call allborn
      if(flg_withdamp) then
         if(valid_emitter(0).or.
     1        valid_emitter(1).or.valid_emitter(2)) then
            call gen_real_phsp_isr(xrad,
     1           jac_over_csi,jac_over_csi_p,jac_over_csi_m,jac_over_csi_s)
            xjac=jac_over_csi*kn_csitilde*kn_csimax**2*kn_jacborn*ww*hc2
            ! xjac=jac_over_csi*kn_csi*kn_csimax*kn_jacborn*ww*hc2
         endif

         rrr_damp_rem_tot=0
         do alr=1,flst_nalr
            rrr_damp_rem_arr(alr)=0
            rrr_rem_arr(alr)=0
         enddo
         do kn_emitter=0,nlegborn
            if(valid_emitter(kn_emitter)) then
               if(kn_emitter.le.2) then
c     No need to generate phase space; it is already available
                  call setscalesbtlreal
c     signed_sig is the sum of all contributions, not including any suppression factor.
c     it is used here only for the NLO analysis.
c     abs_sig is the sum of the absolute value of all contributions,
c     including the suppression factor.
                  call sigreal_damp_rem(xjac,signed_sig,abs_sig,
     &                 rrr_damp_rem_arr,rrr_rem_arr)
                  if(flg_nlotest) then
                     if(flg_detailedNLO) then
                        do alr=1,flst_nalr
                           tmp = rrr_rem_arr(alr)
                           if(pwhg_isfinite(tmp) .and. tmp /= 0) then
c     tmp is zero unless flst_emitter(alr)=kn_emitter
                              flst_currentalr = alr
                              call analysis_driver(tmp,1)
                           endif
                        enddo
                     else
                        if(pwhg_isfinite(signed_sig)
     1                       .and. signed_sig /= 0) then
                           call analysis_driver(signed_sig,1)
                        endif
                     endif
                  endif
                  rrr_damp_rem_tot=rrr_damp_rem_tot+abs_sig
               else
                  call gen_real_phsp_fsr(xrad,
     1                 jac_over_csi,jac_over_csi_coll,jac_over_csi_s)
                  xjac=jac_over_csi*kn_csi*kn_csimax
     2                *kn_jacborn*ww*hc2

                  call setscalesbtlreal
c     signed_sig is the sum of all contributions, not including any suppression factor.
c     it is used here only for the NLO analysis.
c     abs_sig is the sum of the absolute value of all contributions,
c     including the suppression factor.
                  call sigreal_damp_rem(xjac,signed_sig,abs_sig,
     &                 rrr_damp_rem_arr,rrr_rem_arr)
                  if(flg_nlotest) then
                     if(flg_detailedNLO) then
                        do alr=1,flst_nalr
                           tmp = rrr_rem_arr(alr)
                           if(pwhg_isfinite(tmp) .and. tmp /= 0) then
c     tmp is zero unless flst_emitter(alr)=kn_emitter
                              flst_currentalr = alr
                              call analysis_driver(tmp,1)
                           endif
                        enddo
                     else
                        if(pwhg_isfinite(signed_sig)
     1                       .and. signed_sig /= 0) then
                           call analysis_driver(signed_sig,1)
                        endif
                     endif
                  endif
                  rrr_damp_rem_tot=rrr_damp_rem_tot+abs_sig
               endif
            endif
         enddo
      else
         if(flg_debug) then
            write(*,*) 'Sigremnant: Should NEVER get here!'
         endif
         rrr_damp_rem_tot=0
      endif

      if (.not. pwhg_isfinite(rrr_damp_rem_tot)) then
         call increasecnt("NaN's in sigremnant")
         rrr_damp_rem_tot=0
         rrr_damp_rem_arr=0
      endif

      retval = rrr_damp_rem_tot

      end

      subroutine sigreal_damp_rem(xjac,signed_sig,abs_sig,r1,r0)
c     signed_sig is the sum of all contributions, not including any suppression factor.
c     it is used above for the NLO analysis.
c     abs_sig is the sum of the absolute value of all contributions,
c     including the suppression factor.
c     r1 is the array of results for each alr, including the suppression factor, and
c     with their correct sign
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 xjac,suppfact,signed_sig,abs_sig
      real * 8 r0(maxalr),r1(maxalr)
      integer alr
      signed_sig=0
      abs_sig=0
      if (xjac == 0d0) then
         r0 = 0d0
      else
         call sigreal_btl0(r0,1)
      endif
c     now r0 contains non-zero entries only for the alr that have
c     flst_emitter(alr) = kn_emitter
c     In r1 are accumulated all alr.
      do alr=1,flst_nalr
         if(kn_emitter.eq.flst_emitter(alr)) then
            if(kn_emitter.le.2) then
               r0(alr)=r0(alr)/((1-kn_y**2)*(kn_csitilde*kn_csimax)**2)
            else
               r0(alr)=r0(alr)/((1-kn_y)*(kn_csitilde*kn_csimax)**2)
            endif
            r0(alr)=r0(alr)*xjac
            flst_currentalr = alr
            call rmn_suppression(suppfact)
c     The integrated value r1 is suppressed
c     Below saying r1(alr)=r0(alr)*suppfact has the same effect
c     r1 was zeroed before the kn_emitter loop
            r1(alr)=r1(alr)+r0(alr)*suppfact
c     signed_sig is passed to the analysis (no suppression)
            signed_sig=signed_sig+r0(alr)
            abs_sig=abs_sig+abs(r0(alr))*suppfact
         endif
      enddo
      end

      
