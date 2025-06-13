      subroutine btildecoll(xrad,rescoll,www)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      include 'pwhg_par.h'
      real * 8 xrad,rescoll(flst_nborn),www,rescfac
      real * 8 un
      parameter (un=1d0)
c pdfb: born pdf's, pdfs: pdfs with scaled x->x/z
      real * 8
     1     pdfb1(-pdf_nparton:pdf_nparton),
     2     pdfb2(-pdf_nparton:pdf_nparton),
     3     pdfs1(-pdf_nparton:pdf_nparton),
     4     pdfs2(-pdf_nparton:pdf_nparton),
     5     pdfs1sng,pdfs2sng
      real * 8 z,omzpgg,omzpqq,ppqq,omzpqg,ppqg,omzpgq,ppgq,
     1         sb,tot,plfrc,plfr0,z1,z2,xjac1,xjac2,rm1,rm2,res1,res2,
     2         x,xjac
      external omzpgg,omzpqq,ppqq,omzpqg,ppqg,omzpgq,ppgq
      integer j,jb,fl1,fl2
      logical pwhg_isfinite
      real * 8 chargeofparticle,resgam1,resgam2,
     1         ch1,ch2
      external pwhg_isfinite,chargeofparticle
      real * 8 pdfs1sng_gamma,pdfs2sng_gamma
      real * 8 pdfsmb1sng_gamma, pdfsmb2sng_gamma
      real * 8 tmp
c  ---------------------------------------------------------------------
cc ...IB: moved to external functions      
c
cc omz:(1-z)*pgg(z), 2.106 of FNO2007
c      omzpgg(z)=2*ca*(z+(1-z)**2/z+z*(1-z)**2)
cc 2.103
c      omzpqq(z)=cf*(1+z**2)
cc ppqq: (p prime qq) d pqq'/ d epsilon, from 2.103
c      ppqq(z)=-cf*(1-z)
cc 2.104
c      omzpqg(z)=cf*(1+(1-z)**2)/z*(1-z)
c      ppqg(z)=-cf*z
cc 2.105
c      omzpgq(z)=tf*(1-2*z*(1-z))*(1-z)
c      ppgq(z)=-tf*2*z*(1-z)
cc
c-----------------------------------------------------------------------     
c ...IB: similar log coefficients for polarized and unpolarized
c Statement Functions
c log coeff., from 2.102, with partonic s=sb/z
      plfrc(z)=1/(1-z)*log(sb/z/st_mufact2)+2*log(1-z)/(1-z)
c same, with soft limit s
      plfr0(z)=1/(1-z)*log(sb/st_mufact2)+2*log(1-z)/(1-z)    
c-----------------------------------------------------------------------
c End Statement Functions
c If there is no radiation from production, no remnants!
      do j=1,flst_nreson
         if(flst_reslist(j).eq.0) exit
      enddo
      if(j.eq.flst_nreson+1) then
c no radiation from production!
         rescoll=0
         return
      endif

      if(flg_collremnsamp) then
         xjac=630*(1-xrad)**4*xrad**4
         x=xrad**5*(70*xrad**4-315*xrad**3+540*xrad**2-420*xrad+126)
      else
         xjac=1
         x=xrad
      endif
c      z1=1-par_isrtinycsi-(1-kn_xb1-par_isrtinycsi)*x
      z1=1-(1-kn_xb1)*par_isrtinycsi-(1-kn_xb1)*(1-par_isrtinycsi)*x
      xjac1=(1-kn_xb1)*xjac
c      z2=1-par_isrtinycsi-(1-kn_xb2-par_isrtinycsi)*x
      z2=1-(1-kn_xb2)*par_isrtinycsi-(1-kn_xb2)*(1-par_isrtinycsi)*x
      xjac2=(1-kn_xb2)*xjac

      sb=kn_sborn
      if(.not.flg_minlo) then
         rescfac = 1
      endif
c loop over all Born configurations
      tot=0
      do jb=1,flst_nborn
c the following is to be done once if minlo is not set
         if(flg_minlo.or.jb.eq.1) then
c the following also sets st_mufact2 according to the underlying Born index jb,
c sets rescfac to the rescaling factors needed by the virtual, real and collinear remnants
c (if the flag argument 2 is given)
            if(flg_minlo) call setlocalscales(jb,2,rescfac)
c See 7.224 and 7.225 of FNO2007
c Remnant: 1/(1-z)_+=1/(1-z)_(1-x)+log(1-x) delta(1-z)
c Remnant: log(1-z)/(1-z)_+ = log(1-z)/(1-z)_(1-x)+
c          log(1-x)^2/2 delta(1-z)
            if (kn_xb1 < 1) then
               rm1=log(1-kn_xb1)*log(sb/st_mufact2)
     1             + log(1-kn_xb1)**2
            else
               rm1 = 0d0
            endif
            rm2=log(1-kn_xb2)*log(sb/st_mufact2)
     1           + log(1-kn_xb2)**2
c     get pdfs at underlying born x values
            call pdfcall(1,kn_xb1,pdfb1)
            call pdfcall(2,kn_xb2,pdfb2)
c     get pdfs at underlying born x/z value values
            call pdfcall(1,kn_xb1/z1,pdfs1)
            call pdfcall(2,kn_xb2/z2,pdfs2) 
c     Compute the singlet
            pdfs1sng=0
            pdfs2sng=0
            do j=1,st_nlight
               pdfs1sng=pdfs1sng+pdfs1(j)+pdfs1(-j)
               pdfs2sng=pdfs2sng+pdfs2(j)+pdfs2(-j)
            enddo
            if(flg_with_em) then
c     Compute the singlet for gamma induced
               pdfs1sng_gamma=0
               pdfs2sng_gamma=0
               pdfsmb1sng_gamma=0
               pdfsmb2sng_gamma=0
               do j=1,st_nlight
                  pdfs1sng_gamma=
     1                 pdfs1sng_gamma+pdfs1(j)*chargeofparticle(j)**2
     2                 +pdfs1(-j)*chargeofparticle(-j)**2
                  pdfs2sng_gamma=pdfs2sng_gamma+pdfs2(j)
     1                 *chargeofparticle(j)**2
     2                 +pdfs2(-j)*chargeofparticle(-j)**2
                  pdfsmb1sng_gamma=pdfsmb1sng_gamma+
     1                 (pdfs1(j)/z1-pdfb1(j))*  chargeofparticle(j)**2 +
     2                 (pdfs1(-j)/z1-pdfb1(-j))*chargeofparticle(-j)**2
                  pdfsmb2sng_gamma=pdfsmb2sng_gamma+
     1                 (pdfs2(j)/z2-pdfb2(j)) * chargeofparticle(j)**2 +
     2                 (pdfs2(-j)/z2-pdfb2(-j))*chargeofparticle(-j)**2
               enddo
            endif
         endif

         res1    = 0d0
         res2    = 0d0
         resgam1 = 0d0
         resgam2 = 0d0

         fl1=flst_born(1,jb)
         fl2=flst_born(2,jb)

         ch1=chargeofparticle(fl1)
         ch2=chargeofparticle(fl2)

         if(fl1.eq.0) then
c FNO2007, 2.102
c gg remnant
c Remember 1/z factor, 4.23 and 7.206 of FNO2007
            res1=        omzpgg(z1)*plfrc(z1)/z1  * pdfs1(fl1)*xjac1 
     #                 - omzpgg(un)*plfr0(z1)     * pdfb1(fl1)*xjac1
     #                 + omzpgg(un)*rm1           * pdfb1(fl1)
c qg remnant
            res1=res1 + (omzpqg(z1)*plfrc(z1)
     #                             - ppqg(z1))/z1 * pdfs1sng*xjac1 
         elseif(flg_with_em.and.fl1.eq.22) then
c qgamma remnant
            resgam1=resgam1 + (omzpqg(z1)*plfrc(z1)
     #                             - ppqg(z1))/z1 * pdfs1sng_gamma*xjac1
     #                              /cf
            if (pdf_dis_photon) then
                resgam1 = resgam1  + 
     #                    ( (9d0 + 5*z1)/4. + 
     #                      (1d0 + z1**2)/(1d0 - z1)*
     #                        (-0.75d0 + log((1d0 - z1)/z1)) )*
     #                   pdfsmb1sng_gamma
            endif

         else
            if(abs(fl1) <= 6) then
c     qq remnant            
               res1=       (omzpqq(z1)*plfrc(z1)
     1              - ppqq(z1))/z1 * pdfs1(fl1)*xjac1 
     2              - omzpqq(un)*plfr0(z1)     * pdfb1(fl1)*xjac1
     3              + omzpqq(un)*rm1           * pdfb1(fl1)
c WEW, qq remnant with photon
               if(flg_with_em) then
                  resgam1 =   ((omzpqq(z1)*plfrc(z1)
     1                 - ppqq(z1))/z1 * pdfs1(fl1)*xjac1 
     2                 - omzpqq(un)*plfr0(z1)     * pdfb1(fl1)*xjac1
     3                 + omzpqq(un)*rm1           * pdfb1(fl1))
     4                 * ch1**2 /cf
c gammaq DIS subtraction
                  if (pdf_dis_photon) then
                     resgam1 = resgam1  - 
     1                    ( (9d0 + 5*z1)/4. + 
     2                    (1d0 + z1**2)/(1d0 - z1)*
     3                    (-0.75d0 + log((1d0 - z1)/z1)) )*
     4                    (pdfs1(fl1)/z1-pdfb1(fl1)) * ch1**2
                  endif
               endif
c gq remnant
               res1=res1 + (omzpgq(z1)*plfrc(z1)
     1              - ppgq(z1))/z1   * pdfs1(0)*xjac1 
               if(flg_with_em .and. pdf_nparton >= 22) then
c     gamma-q remnant; as gluon-quark remnant divided by tf, times
c     the charge square times 3 for colours in case of quarks
                  resgam1 = resgam1  + (omzpgq(z1)*plfrc(z1)
     1                 - ppgq(z1))/z1   * pdfs1(22)*xjac1 /tf *ch1**2 * 3
               endif
c qgamma DIS subtraction 
               if (flg_with_em.and.pdf_dis_photon) then
                  resgam1 = resgam1  - 3d0*
     1                 ( log( (1d0-z1)/z1 )*( z1**2+(1d0-z1)**2 )
     2                 - 8d0*z1**2 + 8d0*z1 - 1d0  )*
     3                 (pdfs1(fl1)/z1-pdfb1(fl1)) * ch1**2
               endif
            elseif(abs(fl1)>=11) then
               if(flg_with_em .and. pdf_nparton >= 22) then
                  resgam1 = resgam1  + (omzpgq(z1)*plfrc(z1)
     1                 - ppgq(z1))/z1   * pdfs1(22)*xjac1 /tf *ch1**2
               endif
            endif
         endif
c Copy the above, replace 1 with 2 with the editor
         if(fl2.eq.0) then
c FNO2007, 2.102
c gg remnant
            res2=        omzpgg(z2)*plfrc(z2)/z2  * pdfs2(fl2)*xjac2 
     1           - omzpgg(un)*plfr0(z2)     * pdfb2(fl2)*xjac2
     2           + omzpgg(un)*rm2           * pdfb2(fl2)
c qg remnant
            res2=res2 + (omzpqg(z2)*plfrc(z2)
     1           - ppqg(z2))/z2 * pdfs2sng*xjac2 
         elseif(flg_with_em.and.fl2.eq.22) then
c qgamma remnant
            resgam2=resgam2 + (omzpqg(z2)*plfrc(z2)
     #                             - ppqg(z2))/z2 * pdfs2sng_gamma*xjac2
     #                              /cf
c gammaq DIS subtraction 
            if (pdf_dis_photon) then
                resgam2 = resgam2  + 
     #                    ( (9d0 + 5*z2)/4. + 
     #                      (1d0 + z2**2)/(1d0 - z2)*
     #                        (-0.75d0 + log((1d0 - z2)/z2)) )*
     #                   pdfsmb2sng_gamma
            endif

         else
            if(abs(fl2) <= 6) then

c qq remnant            
               res2=       (omzpqq(z2)*plfrc(z2)
     1              - ppqq(z2))/z2 * pdfs2(fl2)*xjac2 
     2              - omzpqq(un)*plfr0(z2)     * pdfb2(fl2)*xjac2
     3              + omzpqq(un)*rm2           * pdfb2(fl2)
               if(flg_with_em) then
c     WEW, qq remnant with photon
                  resgam2 =   ((omzpqq(z2)*plfrc(z2)
     1                 - ppqq(z2))/z2 * pdfs2(fl2)*xjac2 
     2                 - omzpqq(un)*plfr0(z2)     * pdfb2(fl2)*xjac2
     3                 + omzpqq(un)*rm2           * pdfb2(fl2))
     4                 * ch2**2 /cf
c     qq DIS subtraction 
                  if (pdf_dis_photon) then
                     resgam2 = resgam2  - 
     1                    ( (9d0 + 5*z2)/4. + 
     2                    (1d0 + z2**2)/(1d0 - z2)*
     3                    (-0.75d0 + log((1d0 - z2)/z2)) )*
     4                    (pdfs2(fl2)/z2-pdfb2(fl2)) * ch2**2
                  endif
               endif
c gq remnant
               res2=res2 + (omzpgq(z2)*plfrc(z2)
     1                           - ppgq(z2))/z2   * pdfs2(0)*xjac2 

               if(flg_with_em .and. pdf_nparton >= 22) then
c     gamma-q remnant; as gluon-quark remnant divided by tf, times
c     the charge square times 3 for colours in case of quarks
                  resgam2 = resgam2  + (omzpgq(z2)*plfrc(z2)
     1                 - ppgq(z2))/z2   * pdfs2(22)*xjac2 /tf *ch2**2 * 3
               endif
c     qgamma DIS subtraction 
               if (flg_with_em.and.pdf_dis_photon) then
                  resgam2 = resgam2  - 
     1                 ( log( (1d0-z2)/z2 )*( z2**2+(1d0-z2)**2 )
     2                 - 8d0*z2**2 + 8d0*z2 - 1d0  )*
     3                 (pdfs2(fl2)/z2-pdfb2(fl2)) * ch2**2
               endif
            elseif(abs(fl2)>=11) then
               if(flg_with_em .and. pdf_nparton >= 22) then
                  resgam2 = resgam2  + (omzpgq(z2)*plfrc(z2)
     1                 - ppgq(z2))/z2   * pdfs2(22)*xjac2 /tf *ch2**2
               endif
            endif
         endif
         if (flg_QEDonly) then
            res1=0
            res2=0
         endif
         rescoll(jb)=( res1*pdfb2(fl2)+res2*pdfb1(fl1) )
     #    *br_born(jb)*st_alpha/(2*pi)*kn_jacborn * rescfac
         if(flg_with_em) then
c     WEW, add photon remnants       
            rescoll(jb)=rescoll(jb)
     1           +(resgam1*pdfb2(fl2)+resgam2*pdfb1(fl1))
     2           *br_born(jb)*em_alpha/(2*pi)*kn_jacborn * rescfac
         endif
         tot=tot+rescoll(jb)
      enddo
      if (.not.pwhg_isfinite(tot)) then
         do jb=1,flst_nborn
            rescoll(jb)=0d0
         enddo
         tot=0d0
      endif
      if(flg_nlotest) then
         if(flg_detailedNLO) then
            do jb=1,flst_nborn
c     contributions not belonging to the current resonance group turn out
c     to be zero anyhow, so the following if isn't really necessary ...
               if(flst_ibornresgroup == flst_bornresgroup(jb)) then
                  flst_currentuborn = jb
                  tmp = rescoll(jb)*www
                  if(pwhg_isfinite(tmp) .and. tmp /= 0) then
                     call analysis_driver(tmp,0)
                  endif
               endif
            enddo
         else
            tot=tot*www
            call analysis_driver(tot,0)
         endif
      endif
      end

C-----------------------------------------------------------------------
C ...IB: define splitting and derivatives as external functtions
c
      function omzpgg(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 omzpgg,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)    
cc omz:(1-z)*pgg(z), 2.106 of FNO2007
      omzpgg = 2*ca*(z+(1-z)**2/z+z*(1-z)**2)*(1-ipol)
     &       + 2*ca*(1+(1-z)-2*z*(1-z))*ipol
      end

      function omzpqq(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 omzpqq,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)   
c 2.103
      omzpqq = cf*(1+z**2) !same for unpolarized and polarized case
      end

      function ppqq(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 ppqq,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)  
c ppqq: (p prime qq) d pqq'/ d epsilon, from 2.103
      ppqq = -cf*(1-z)     !same for unpolarized and polarized case
      end

      function omzpqg(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 omzpqg,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)  
c 2.104
      omzpqg = cf*(1+(1-z)**2)/z*(1-z)*(1-ipol)
     &       + cf*(1-(1-z)**2)/z*(1-z)*ipol
      end

      function ppqg(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 ppqg,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)  
c
      ppqg = -cf*z*(1-ipol)
     &     +  cf*2*(1-z)*ipol
      end

      function omzpgq(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 omzpgq,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)
c     2.105
      omzpgq = tf*(1-2*z*(1-z))*(1-z)*(1-ipol)
     &       + tf*(2*z-1)*(1-z)*ipol
      end

      function ppgq(z)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 ppgq,z
      integer ipol
c     ipol = 1 if flg_polarised_dis = true, 0 otherwise
      ipol = merge(1, 0, flg_polarised_dis)  
c
      ppgq = -tf*2*z*(1-z)*(1-ipol)
     &       -tf*2*(1-z)*ipol
      end
cc
C-----------------------------------------------------------------------

      function chargeofparticle(id)
c WEW,  chargeofparticle() added
c Returns the electric charge (in units of the positron charge)
c of particle id (pdg conventions, gluon is zero)
      implicit none
      real * 8 chargeofparticle
      integer id
      if(abs(id).gt.0.and.abs(id).le.6) then
         if(mod(abs(id),2).eq.0) then
            chargeofparticle =  2d0/3
         else
            chargeofparticle = -1d0/3
         endif
      elseif(abs(id).gt.10.and.abs(id).le.16) then
         if(mod(abs(id),2).ne.0) then
            chargeofparticle = -1d0
         else
            chargeofparticle = 0
         endif
      else
         chargeofparticle=0
      endif
      if(id<0) chargeofparticle = - chargeofparticle
      end
