c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine esnbnd2  --  atom-by-atom non-bonded 12-10  Hessian  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "esnbnd2" calculates second derivatives of non-bonded 12-10 terms
c     for one atom at the time. 
c
c
      subroutine esnbnd2 (iatom)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'vdwpot.i'
      include 'snbnb.i'
      integer iatom,jcell
      integer i,j,k,ii
      integer exn,exm
      real*8 eps,e,vscale,fgrp
      real*8 rv,rvn,rvm
      real*8 xr,yr,zr
      real*8 rik,rikn,rikm
      real*8 rik2,rik3,rik4,rik5
      real*8 rdrv,nm,mnm
      real*8 ppa,dpa,ddpa
      real*8 ppb,dpb,ddpb
      real*8 dpc,ddpc
      real*8 taper,dtaper,d2taper
      real*8 dedr,d2edr2,de,d2e
      real*8 d2edx,d2edy,d2edz
      real*8 term(3,3)
      character*6 mode
      logical proceed
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     check that iatom has 12-10 non-bonded interactions
c
      proceed = .false.
      do i = 1, nsnb
         if (isnb(1,i).eq.0 .and. isnb(2,i).eq.0) exit
         if (isnb(1,i).eq.iatom .or. isnb(2,i).eq.iatom) then
            proceed = .true.
            exit
         endif
      enddo
      if (.not. proceed) return
c
c     find the Hessian elements for involved atoms
c
      do ii = 1, nsnb
         if (isnb(1,ii).eq.0 .and. isnb(2,ii).eq.0) then
            exit
         else if (isnb(1,ii).eq.iatom) then
            i = isnb(1,ii)
            k = isnb(2,ii)
         else if (isnb(2,ii).eq.iatom) then
            i = isnb(2,ii)
            k = isnb(1,ii)
         else
            cycle
         end if
         rv = rsnd(ii)
         eps = epssnb(ii)
         vscale = 1.0d0
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (k .ne. i)
c
c     compute the value of the distance deviation
c
         if (proceed) then
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            call image (xr,yr,zr)
            rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
            if (rik2 .le. off2) then
               rik = sqrt(rik2)
c
c     set interaction scaling coefficients for connected atoms
c
               do j = 1, n12(i)
                  if (i12(j,i) .eq. k) vscale = v2scale
               end do
               do j = 1, n13(i)
                  if (i13(j,i) .eq. k) vscale = v3scale
               end do
               do j = 1, n14(i)
                  if (i14(j,i) .eq. k) vscale = v4scale
               end do
               do j = 1, n15(i)
                  if (i15(j,i) .eq. k) vscale = v5scale
               end do
               eps = eps * vscale
c
c     calculate the contribution to the Hessian
c
               exn = 10
               exm = 2
               rvn = rv**exn
               rvm = rv**exm
               rikn = rik**exn
               rikm = rik**exm
               rdrv = rik + dhal*rv
               ppa= (dhal + 1.0d0) / rdrv
               nm = exn / exm
               mnm = (exm + exn) / exm
               ppb = rikm + ghal*rvm

               e = eps * rvn * ppa**exn
     &             * (nm * (1.0d0+ghal) * rvm/ppb - mnm)

               dpa = e / rdrv
               dpb = eps * rvn * rvm * ppa**exn * rikm/rik
               dpc = (1.0d0 + ghal) / ppb**2

               dedr = -exn * (dpa + dpb*dpc)

               ddpa = dedr/rdrv - e/rdrv**2
               ddpb = dpb * ((exm-1.0d0)/rik - exn/rdrv)
               ddpc = -2.0d0 * exm * rikm * dpc / (rik * ppb)

               d2edr2 = -exn * (ddpa + ddpb*dpc + dpb*ddpc)
c
c     use energy switching if near the cutoff distance
c
               if (rik2 .gt. cut2) then
                  rik3 = rik2 * rik
                  rik4 = rik2 * rik2
                  rik5 = rik3 * rik2
                  taper = c5*rik5 + c4*rik4 + c3*rik3
     &                       + c2*rik2 + c1*rik + c0
                  dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                        + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                  d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                         + 6.0d0*c3*rik + 2.0d0*c2
                  d2edr2 = e*d2taper + 2.0d0*dedr*dtaper + d2edr2*taper
                  dedr = e*dtaper + dedr*taper
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedr = dedr * fgrp
                  d2edr2 = d2edr2 * fgrp
               end if
c
c     get chain rule terms for van der Waals Hessian elements
c
               de = dedr / rik
               d2e = (d2edr2-de) / rik2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nsnb
         if (isnb(1,ii).eq.0 .and. isnb(2,ii).eq.0) then
            exit
         else if (isnb(1,ii).eq.iatom) then
            i = isnb(1,ii)
            k = isnb(2,ii)
         else if (isnb(2,ii).eq.iatom) then
            i = isnb(2,ii)
            k = isnb(1,ii)
         else
            cycle
         end if
         rv = rsnd(ii)
         eps = epssnb(ii)
         vscale = 1.0d0
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the value of the distance deviation
c
         if (proceed) then
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
c
c     set interaction scaling coefficients for connected atoms
c
            do j = 1, n12(i)
               if (i12(j,i) .eq. k) vscale = v2scale
            end do
            do j = 1, n13(i)
               if (i13(j,i) .eq. k) vscale = v3scale
            end do
            do j = 1, n14(i)
               if (i14(j,i) .eq. k) vscale = v4scale
            end do
            do j = 1, n15(i)
               if (i15(j,i) .eq. k) vscale = v5scale
            end do
c
c     loop on cells
c
            do jcell = 1, ncell
               call imager (xr,yr,zr,jcell)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  if (use_polymer) then
                     if (rik2 .le. polycut2) then
                        eps = eps * vscale
                     endif
                  endif
                  exn = 10
                  exm = 2
                  rvn = rv**exn
                  rvm = rv**exm
                  rikn = rik**exn
                  rikm = rik**exm
                  rdrv = rik + dhal*rv
                  ppa= (dhal + 1.0d0) / rdrv
                  nm = exn / exm
                  mnm = (exm + exn) / exm
                  ppb = rikm + ghal*rvm
   
                  e = eps * rvn * ppa**exn
     &                * (nm * (1.0d0+ghal) * rvm/ppb - mnm)

                  dpa = e / rdrv
                  dpb = eps * rvn * rvm * ppa**exn * rikm/rik
                  dpc = (1.0d0 + ghal) / ppb**2
   
                  dedr = -exn * (dpa + dpb*dpc)

                  ddpa = dedr/rdrv - e/rdrv**2
                  ddpb = dpb * ((exm-1.0d0)/rik - exn/rdrv)
                  ddpc = -2.0d0 * exm * rikm * dpc / (rik * ppb)

                  d2edr2 = -exn * (ddpa + ddpb*dpc + dpb*ddpc)
c
c     use energy switching if near the cutoff distance
c
               if (rik2 .gt. cut2) then
                  rik3 = rik2 * rik
                  rik4 = rik2 * rik2
                  rik5 = rik3 * rik2
                  taper = c5*rik5 + c4*rik4 + c3*rik3
     &                       + c2*rik2 + c1*rik + c0
                  dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                        + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                  d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                         + 6.0d0*c3*rik + 2.0d0*c2
                  d2edr2 = e*d2taper + 2.0d0*dedr*dtaper + d2edr2*taper
                  dedr = e*dtaper + dedr*taper
               end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     dedr = dedr * fgrp
                     d2edr2 = d2edr2 * fgrp
                  end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                  de = dedr / rik
                  d2e = (d2edr2-de) / rik2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
         end if
      end do
      return
      end subroutine
