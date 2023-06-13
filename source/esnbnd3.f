c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine esnbnd3  -- calculates the energy of 12-10 non-bonded terms    ##
c     ##                                                                            ##
c     ################################################################################
c
c
c     "esnbnd" calculates the energy for all non-bonded 12-10 terms and 
c     partitions the energy among the atoms.
c
c
      subroutine esnbnd3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
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
      real*8 taper,fgrp
      real*8 xr,yr,zr
      character*6 mode
      logical proceed
      logical header,huge
c
c
c     zero out the energy and partitioning terms
c
      nesnb = 0
      esnb = 0.0d0
      do i = 1, n
         aesnb(i) = 0.0d0
      end do
      header = .true.
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
               ppa= (dhal + 1.0d0) / (rik + dhal*rv)
               nm = exn / exm
               mnm = (exm + exn) / exm
               ppb = rikm + ghal*rvm
               e = eps * rvn * ppa**exn
     &             * (nm * (1.0d0+ghal) * rvm/ppb - mnm)
c
c     use energy switching if near the cutoff distance
c
               if (rik2 .gt. cut2) then
                  rik3 = rik2 * rik
                  rik4 = rik2 * rik2
                  rik5 = rik2 * rik3
                  taper = c5*rik5 + c4*rik4 + c3*rik3
     &                       + c2*rik2 + c1*rik + c0
                  e = e * taper
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the overall 12-10 energy components
c
               if (e .ne. 0.0d0) then
                  nesnb = nesnb + 1
                  esnb = esnb + e
                  aesnb(i) = aesnb(i) + 0.5d0*e
                  aesnb(k) = aesnb(k) + 0.5d0*e
               end if
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 10.0d0)
               if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual vdW 12-10 Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          20x,'Minimum',4x,'Actual',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20) i,name(i),k,name(k),
     &                                  rv,rik,e
   20             format (' 12-10vdW',3x,2(i7,'-',a3),
     &                    13x,2f10.4,f12.4)
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
                  ppa= (dhal + 1.0d0) / (rik + dhal*rv)
                  nm = exn / exm
                  mnm = (exm + exn) / exm
                  ppb = rikm + ghal*rvm
                  e = eps * rvn * ppa**exn
     &               * (nm * (1.0d0+ghal) * rvm/ppb - mnm)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall 12-10 energy components
c
                  if (e .ne. 0.0d0) then
                     nesnb = nesnb + 1
                     if (i .eq. k) then
                        esnb = esnb + 0.5d0*e
                        aesnb(i) = aesnb(i) + 0.5d0*e
                     else
                        esnb = esnb + e
                        aesnb(k) = aesnb(k) + 0.5d0*e
                        aesnb(k) = aesnb(k) + 0.5d0*e
                     end if
                  end if
c
c     increment the total intermolecular energy
c
                  einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
                     end if
                     write (iout,40) i,name(i),k,name(k),
     &                               rv,rik,e
   40                format (' 12-12vdW',3x,2(i7,'-',a3),3x,
     &                       '(XTAL)',4x,2f10.4,f12.4)
                  end if
               end if
            end do
         end if
      end do
      return
      end subroutine
