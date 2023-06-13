c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine esnbnd1  -- energy of 12-10 non-bonded terms & derivarives   ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "esnbnd" calculates the energy and first derivatives with respect to
c     Cartesian coordinates for all non-bonded 12-10 terms.
c
c
      subroutine esnbnd1
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdwpot.i'
      include 'snbnb.i'
      integer i,j,k,ii
      integer exn,exm
      real*8 e,eps,vscale
      real*8 rv,rvn,rvm
      real*8 rik,rikn,rikm
      real*8 rik2,rik3
      real*8 rik4,rik5
      real*8 ppa,ppb,nm,mnm
      real*8 taper,fgrp,rdrv
      real*8 xr,yr,zr
      real*8 dpa,dpb,dpc
      real*8 de,dedr,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      character*6 mode
      logical proceed
c
c
c     zero out the energy and derivatives
c
      esnb = 0.0d0
      do i = 1, n
         desnb(1,i) = 0.0d0
         desnb(2,i) = 0.0d0
         desnb(3,i) = 0.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     compute evergy the 12-10 non-bonded energy term
c
      do ii = 1, nsnb
         i = isnb(1,ii)
         k = isnb(2,ii)
         rv = rsnd(ii)
         eps = epssnb(ii)
         vscale = 1.0d0
c
c     decide whether to compute the current interaction
c
         if (i .eq. k) cycle
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (use(i) .or. use(k))
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
c     get the interaction energy
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
c
c     use energy switching if near the cutoff distance
c
               if (rik2 .gt. cut2) then
                  rik3 = rik2 * rik
                  rik4 = rik2 * rik2
                  rik5 = rik2 * rik3
                  taper = c5*rik5 + c4*rik4 + c3*rik3
     &                       + c2*rik2 + c1*rik + c0
                  dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                        + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                  dedr = e*dtaper + dedr*taper
                  e = e * taper
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  then
                  e = e * fgrp
                  dedr = dedr * fgrp
               end if
c
c     find the chain rule terms for derivative components
c
               de = dedr / rik
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the overall 12-10 energy components
c
               esnb = esnb + e
               desnb(1,i) = desnb(1,i) + dedx
               desnb(2,i) = desnb(2,i) + dedy
               desnb(3,i) = desnb(3,i) + dedz
               desnb(1,k) = desnb(1,k) - dedx
               desnb(2,k) = desnb(2,k) - dedy
               desnb(3,k) = desnb(3,k) - dedz
c
c     increment the internal virial tensor components
c
               vxx = xr * dedx
               vyx = yr * dedx
               vzx = zr * dedx
               vyy = yr * dedy
               vzy = zr * dedy
               vzz = zr * dedz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
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
         i = isnb(1,ii)
         k = isnb(2,ii)
         rv = rsnd(ii)
         eps = epssnb(ii)
         vscale = 1.0d0
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (use(i) .or. use(k))
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
            do j = 1, ncell
               call imager (xr,yr,zr,j) 
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
c
c     get the interaction energy
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
     &               * (nm * (1.0d0+ghal) * rvm/ppb - mnm)

                  dpa = e / rdrv
                  dpb = eps * rvn * rvm * ppa**exn * rikm/rik
                  dpc = (1.0d0 + ghal) / ppb**2

                  dedr = -exn * (dpa + dpb*dpc)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                          + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     dedr = e*dtaper + dedr*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  then
                     e = e * fgrp
                     dedr = dedr * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                   de = dedr / rik
                   dedx = de * xr
                   dedy = de * yr
                   dedz = de * zr
c
c     increment the overall 12-10 energy components
c
                   if (i .eq. k)  e = 0.5d0 * e
                   esnb = esnb + e
                   desnb(1,i) = desnb(1,i) + dedx
                   desnb(2,i) = desnb(2,i) + dedy
                   desnb(3,i) = desnb(3,i) + dedz
                   if (i .eq. k) then
                      desnb(1,k) = desnb(1,k) - dedx
                      desnb(2,k) = desnb(2,k) - dedy
                      desnb(3,k) = desnb(3,k) - dedz
                   end if
c
c     increment the internal virial tensor components
c
                   vxx = xr * dedx
                   vyx = yr * dedx
                   vzx = zr * dedx
                   vyy = yr * dedy
                   vzy = zr * dedy
                   vzz = zr * dedz
                   vir(1,1) = vir(1,1) + vxx
                   vir(2,1) = vir(2,1) + vyx
                   vir(3,1) = vir(3,1) + vzx
                   vir(1,2) = vir(1,2) + vyx
                   vir(2,2) = vir(2,2) + vyy
                   vir(3,2) = vir(3,2) + vzy
                   vir(1,3) = vir(1,3) + vzx
                   vir(2,3) = vir(2,3) + vzy
                   vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                   einter = einter + e
               end if
            end do
         end if
      end do
      return
      end subroutine
