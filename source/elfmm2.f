c
c
c     ##################################################################
c     ##  COPYRIGHT (C)  2014  by Marco Foscato & Jay William Ponder  ##
c     ##              All Rights Reserved                             ##
c     ##################################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine elfmm2  -- atom-by-atom LFMM Hessian   ##
c     ##                                                    ##
c     ########################################################
c
c
c     "elfmm2" calls the routines that calculate second derivatives 
c     of all LFMM terms for one atom at the time. 
c
c
      subroutine elfmm2 (i)
      implicit none
      integer i
c
c     Calculate Metal-Ligand Stretch term
c
      call harm4lfmm2 (i)
      call morse4lfmm2 (i)
c
c     Calculate Ligand-Ligand Repulsion term
c
      call vdw4lfmm2 (i)
      call ll4lfmm2 (i)
c
c     Calculate Electron Pairing energy
c
      call pair4lfmm2 (i)
c
c     Calculate Ligand Filed Stabilization term
c
      call lfse4lfmm2 (i)

      return
      end subroutine
c
c
c     ############################################################################
c     ##                                                                        ##
c     ##  subroutine harm4lfmm2  --  atom-by-atom harmonic M-L stretch Hessian  ##   
c     ##                                                                        ##
c     ############################################################################
c
c
c     "harm4lfmm2" calculates second derivative of the LFMM metal-ligand
c     harmonic bond stretching energy for a single atom at a time.
c
c
      subroutine harm4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'hessn.i'
      include 'lfmm.i'
      integer i,j,k,ia,ib
      real*8 ideal,force
      real*8 de,deddt,d2eddt2
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 dt,term
      real*8 termx,termy,termz
      real*8 d2e(3,3)
      logical proceed
c
c
c     compute the Hessian elements 
c
      do k = 1, nhalf
         proceed = .false.
         if (ihalf(1,k) .eq. i) then
            ia = ihalf(1,k)
            ib = ihalf(2,k)
            proceed = .true.
         else if (ihalf(2,k) .eq. i) then
            ia = ihalf(2,k)
            ib = ihalf(1,k)
            proceed = .true.
         end if
         if (proceed) then
            ideal = lhalf(k)
            force = khalf(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
            if (proceed) then
               xab = x(ia) - x(ib)
               yab = y(ia) - y(ib)
               zab = z(ia) - z(ib)
               rab2 = xab*xab + yab*yab + zab*zab
               rab = sqrt(rab2)
               dt = rab - ideal
c
c     Harmonic interaction has functional form U(r) = 1/2 * k * (r - r0)**2
c
               deddt = force*dt
               d2eddt2 = force
c
c     compute chain rule terms needed for derivatives
c
               if (rab .eq. 0.0d0) then
                  de = 0.0d0
                  term = 0.0d0
               else
                  de = deddt / rab
                  term = (d2eddt2-de) / rab2
               end if
               termx = term * xab
               termy = term * yab
               termz = term * zab
               d2e(1,1) = termx*xab + de
               d2e(1,2) = termx*yab
               d2e(1,3) = termx*zab
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yab + de
               d2e(2,3) = termy*zab
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,ia) = hessx(j,ia) + d2e(1,j)
                  hessy(j,ia) = hessy(j,ia) + d2e(2,j)
                  hessz(j,ia) = hessz(j,ia) + d2e(3,j)
                  hessx(j,ib) = hessx(j,ib) - d2e(1,j)
                  hessy(j,ib) = hessy(j,ib) - d2e(2,j)
                  hessz(j,ib) = hessz(j,ib) - d2e(3,j)
               end do
            end if
         endif
      end do
      return
      end subroutine
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine morse4lfmm2  --  atom-by-atom Morse M-L stretch Hessian  ##   
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "morse4lfmm2" calculates second derivative of the LFMM metal-ligand
c     Morse bond stretching energy for a single atom at a time.
c
c
      subroutine morse4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'hessn.i'
      include 'lfmm.i'
      integer i,j,k,ia,ib
      real*8 ideal,bde,alpha
      real*8 expterm
      real*8 dt,deddt,d2eddt2
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 lw,up,lwc,upc
      real*8 lwc2,lwc3,upc2,upc3
      real*8 de,term
      real*8 termx,termy,termz
      real*8 d2e(3,3)
      logical proceed
c
c     compute the Hessian elements of the bond stretching
c
      do k = 1, nmolf
         proceed = .false.
         if (imolf(1,k) .eq. i) then
            ia = imolf(1,k)
            ib = imolf(2,k)
            proceed = .true.
         else if (imolf(2,k) .eq. i) then
            ia = imolf(2,k)
            ib = imolf(1,k)
            proceed = .true.
         end if
         if (proceed) then
            ideal = lmolf(k)
            bde = bmolf(k)
            alpha = amolf(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
            if (proceed) then
               xab = x(ia) - x(ib)
               yab = y(ia) - y(ib)
               zab = z(ia) - z(ib)
               rab2 = xab*xab + yab*yab + zab*zab
               rab = sqrt(rab2)
               dt = rab - ideal
c
c     Morse potential uses energy = BDE * ((1 - e**(-alpha*dt))**2 - 1)
c
               expterm = exp(-alpha*dt)
               deddt = 2.0d0*bde*alpha*(1.0d0-expterm)*expterm
               d2eddt2 = -2.0d0*bde*alpha**2.0d0*expterm
     &                    * (1.0d0-2.0d0*expterm)
c
c     add restrains
c
               lw = ideal - lowmolf
               up = ideal + upmolf
               lwc = lw*lw - rab2
               upc = rab2 - up*up
               upc2 = upc*upc
               upc3 = upc2*upc
               lwc2 = lwc*lwc
               lwc3 = lwc2*lwc
               if (lwc3 .gt. 0.0d0) then
                  deddt = deddt - 6.0d0*w8tmolf*rab*lwc2
                  d2eddt2 = d2eddt2 + 24.0d0*w8tmolf*rab2*lwc 
     &                      - 6.0d0*w8tmolf*lwc2
               endif
               if (upc3 .gt. 0.0d0) then
                  deddt = deddt + 6.0d0*w8tmolf*rab*upc2
                  d2eddt2 = d2eddt2 + 24.0d0*w8tmolf*rab2*upc
     &                      + 6.0d0*w8tmolf*upc2
               endif
c
c     set the chain rule terms for the Hessian elements
c
               if (rab2 .eq. 0.0d0) then
                  de = 0.0d0
                  term = 0.0d0
               else
                  de = deddt / rab
                  term = (d2eddt2-de) / rab2
               end if
               termx = term * xab
               termy = term * yab
               termz = term * zab
               d2e(1,1) = termx*xab + de
               d2e(1,2) = termx*yab
               d2e(1,3) = termx*zab
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yab + de
               d2e(2,3) = termy*zab
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,ia) = hessx(j,ia) + d2e(1,j)
                  hessy(j,ia) = hessy(j,ia) + d2e(2,j)
                  hessz(j,ia) = hessz(j,ia) + d2e(3,j)
                  hessx(j,ib) = hessx(j,ib) - d2e(1,j)
                  hessy(j,ib) = hessy(j,ib) - d2e(2,j)
                  hessz(j,ib) = hessz(j,ib) - d2e(3,j)
               end do
            endif
         endif
      enddo
      return
      end subroutine
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine vdw4lfmm2  --  atom-by-atom VDW L-L interaction Hessian  ##   
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "vdw4lfmm2" calculates second derivative of the LFMM ligand-ligand
c     van der Waals-like interaction for a single atom at a time.
c
c
      subroutine vdw4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'lfmm.i'
      include 'hessn.i'
      integer i,j,k,ia,ib,m
      real*8 a,b
      real*8 dedr,d2edr2
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 de,term
      real*8 termx,termy,termz
      real*8 d2e(3,3)
      logical proceed
c
c
c     compute the Hessian elements 
c
      do k = 1, nvwlf
         proceed = .false.
         if (ivwlf(1,k) .eq. i) then
            ia = ivwlf(1,k)
            ib = ivwlf(2,k)
            proceed = .true.
         else if (ivwlf(2,k) .eq. i) then
            ia = ivwlf(2,k)
            ib = ivwlf(1,k)
            proceed = .true.
         end if
         if (proceed) then
           m = mvwlf(k)
           a = avwlf(k)
           b = bvwlf(k)
c
c     decide whether to compute the current interaction
c
           proceed = .true.
           proceed = (use(ia) .or. use(ib))
c
c     compute the value of the L-L length 
c
           if (proceed) then
              xab = x(ia) - x(ib)
              yab = y(ia) - y(ib)
              zab = z(ia) - z(ib)
              rab2 = xab*xab + yab*yab + zab*zab
              rab = sqrt(rab2)
c
c     Ligand-Ligand interaction as Lennard-Jones-like term
c     U(r) = a/r**m-b/r**6
c
               dedr = -a*m*rab**(-m-1.0d0) + 6.0d0*b*(1.0d0/rab**7.0d0)
               d2edr2 = a*m*(m+1.0d0)*rab**(-m-2.0d0) 
     &                  - 42.0d0*b*rab**(-8.0d0)
c
c     compute chain rule terms needed for derivatives
c
               if (rab .eq. 0.0d0) then
                  de = 0.0d0
                  term = 0.0d0
               else
                  de = dedr / rab
                  term = (d2edr2-de) / rab2
               end if
               termx = term * xab
               termy = term * yab
               termz = term * zab
               d2e(1,1) = termx*xab + de
               d2e(1,2) = termx*yab
               d2e(1,3) = termx*zab
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yab + de
               d2e(2,3) = termy*zab
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,ia) = hessx(j,ia) + d2e(1,j)
                  hessy(j,ia) = hessy(j,ia) + d2e(2,j)
                  hessz(j,ia) = hessz(j,ia) + d2e(3,j)
                  hessx(j,ib) = hessx(j,ib) - d2e(1,j)
                  hessy(j,ib) = hessy(j,ib) - d2e(2,j)
                  hessz(j,ib) = hessz(j,ib) - d2e(3,j)
               end do
            end if
         endif
      end do
      return
      end subroutine
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine ll4lfmm2  --  atom-by-atom L-L pure repulsion Hessian  ##   
c     ##                                                                    ##
c     ########################################################################
c
c
c     "ll4lfmm2" calculates second derivative of the LFMM ligand-ligand
c     pure repulsive energy for a single atom at a time.
c
c
      subroutine ll4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'hessn.i'
      include 'lfmm.i'
      integer i,j,k,ia,ib
      real*8 aaa,enne
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 de,dedr,d2edr2
      real*8 term
      real*8 termx,termy,termz
      real*8 d2e(3,3)
      logical proceed
c
c
c     compute the Hessian elements 
c
      do k = 1, nlllf
         proceed = .false.
         if (illlf(1,k) .eq. i) then
            ia = illlf(1,k)
            ib = illlf(2,k)
            proceed = .true.
         else if (illlf(2,k) .eq. i) then
            ia = illlf(2,k)
            ib = illlf(1,k)
            proceed = .true.
         end if
         if (proceed) then
            enne = mlllf(k)
            aaa = alllf(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            proceed = (use(ia) .or. use(ib))
c
c     compute the value of the L-L distance
c
            if (proceed) then
               xab = x(ia) - x(ib)
               yab = y(ia) - y(ib)
               zab = z(ia) - z(ib)
               rab2 = xab*xab + yab*yab + zab*zab
               rab = sqrt(rab2)
c
c     Ligand-Ligand interaction has functional form U(r) = A*r**(-n) .
c
               dedr = -enne*aaa*rab**(-enne-1)
               d2edr2 = aaa*enne*(enne+1)*rab**(-enne-2)
c
c     compute chain rule terms needed for derivatives
c
               if (rab .eq. 0.0d0) then
                  de = 0.0d0
                  term = 0.0d0
               else
                  de = dedr / rab
                  term = (d2edr2-de) / rab2
               end if
               termx = term * xab
               termy = term * yab
               termz = term * zab
               d2e(1,1) = termx*xab + de
               d2e(1,2) = termx*yab
               d2e(1,3) = termx*zab
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yab + de
               d2e(2,3) = termy*zab
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,ia) = hessx(j,ia) + d2e(1,j)
                  hessy(j,ia) = hessy(j,ia) + d2e(2,j)
                  hessz(j,ia) = hessz(j,ia) + d2e(3,j)
                  hessx(j,ib) = hessx(j,ib) - d2e(1,j)
                  hessy(j,ib) = hessy(j,ib) - d2e(2,j)
                  hessz(j,ib) = hessz(j,ib) - d2e(3,j)
               end do
            end if
         endif
      end do
      return
      end subroutine
c
c
c     ###########################################################################
c     ##                                                                       ##
c     ##  subroutine pair4lfmm2  --  atom-by-atom empiric pairing el. Hessian  ##   
c     ##                                                                       ##
c     ###########################################################################
c
c
c     "pair4lfmm2" calculates second derivative of the LFMM empiric electron
c     pairing energy, computed on the basis of metal-ligand bonds,
c     for a single atom at a time.
c
c
      subroutine pair4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'lfmm.i'
      include 'hessn.i'
      integer i,j,k,ia,ib
      real*8 xab,yab,zab
      real*8 w,rab,rab2
      real*8 de,dedr,d2edr2
      real*8 term
      real*8 termx,termy,termz
      real*8 d2e(3,3)
      real*8 eprm(7)
      logical proceed
c
c
c     compute the Hessian elements 
c
      do k = 1, npair
         proceed = .false.
         if (ipair(1,k) .eq. i) then
            ia = ipair(1,k)
            ib = ipair(2,k)
            proceed = .true.
         else if (ipair(2,k) .eq. i) then
            ia = ipair(2,k)
            ib = ipair(1,k)
            proceed = .true.
         end if
         if (proceed) then
            eprm = ppair(k,:)
            w = wpair(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            proceed = (use(ia) .or. use(ib))
c
c     compute the value of M-L bond length
c
            if (proceed) then
               xab = x(ia) - x(ib)
               yab = y(ia) - y(ib)
               zab = z(ia) - z(ib)
               rab2 = xab*xab + yab*yab + zab*zab
               rab = sqrt(rab2)
c
c     derivatives of polynomial
c
               dedr = w*eprm(2)
               d2edr2 = 0.0d0
               do j = 3, 7
                  dedr = dedr - w*(j-1)*eprm(j)*(rab**(-j))
                  d2edr2 = d2edr2 + w*(j-1)*j*eprm(j)*(rab**(-j-1))
               enddo
c
c     compute chain rule terms needed for derivatives
c
               if (rab .eq. 0.0d0) then
                  de = 0.0d0
                  term = 0.0d0
               else
                  de = dedr / rab
                  term = (d2edr2-de) / rab2
               end if
               termx = term * xab
               termy = term * yab
               termz = term * zab
               d2e(1,1) = termx*xab + de
               d2e(1,2) = termx*yab
               d2e(1,3) = termx*zab
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yab + de
               d2e(2,3) = termy*zab
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,ia) = hessx(j,ia) + d2e(1,j)
                  hessy(j,ia) = hessy(j,ia) + d2e(2,j)
                  hessz(j,ia) = hessz(j,ia) + d2e(3,j)
                  hessx(j,ib) = hessx(j,ib) - d2e(1,j)
                  hessy(j,ib) = hessy(j,ib) - d2e(2,j)
                  hessz(j,ib) = hessz(j,ib) - d2e(3,j)
               end do
            end if
         endif
      end do
      return
      end subroutine
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine lfse4lfmm2  --  atom-by-atom ligand field Hessian  ##   
c     ##                                                                ##
c     ####################################################################
c
c
c     "lfse4lfmm2" is the entry point for the calculation of the
c     second derivatives of ligand field stabilization energy 
c     (according to LFMM method)
c     with respect to Cartesian coordinates by means of finite 
c     difference method. 
c
c     The switch 'dosdfdout'
c     allows to choose between using the finite difference method 
c     available within ICCG routines or having Tinker do the job
c     and collecting only the gradient from ICCG routines.
c     This last possibility is meant for numerical comparison
c     or detailed tests and therefore nearly always ignored.
c
c
      subroutine lfse4lfmm2 (i)
      implicit none
      include 'sizes.i'
      include 'lfmmset.i'
      integer i

      if (dosdfdout) then
         call lfse4lfmm2b (i)
      else
         call lfse4lfmm2a (i)      
      endif
      return
      end subroutine
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine lfse4lfmm2a  --  fast atom-by-atom ligand field Hessian  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "lfse4lfmm2a" obtains Hessian entries. At first call, it calculates
c     and stores all the lfse contributions to the Hessian; the
c     next calls only extract the proper value from the stored data
c     and increment the proper entry of the Hessian. 
c
c
      subroutine lfse4lfmm2a (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'deriv.i'
      include 'hessn.i'
      include 'lfsehess.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      include 'iounit.i'
      integer i,j,k,l,m
      integer ia,ib,ca,cb 
      integer size
      character*4 pia,pib,pca,pcb
      character*16 pt,pt2,blank16
      logical found,proceed,calchess
c
c
c     Identify all metal centers related with atom i
c
      blank16 = '                '
      size = 4
c
c     update list of up-to-date lfse contribs to Hessian
c
      if (updtlfshess) then
         do k = 1, maxlfmm
            sdlfsedone(k) = .false.
         enddo
         updtlfshess = .false.
      endif
c
c        skip d0 and d10
c
      do k = 1, nlfse
         if (nelelfse(k).eq.0 .or. nelelfse(k).eq.10) then
            cycle
         endif
c
c        decide whether to compute the current term
c
         proceed = .false.
         calchess = .false.
         do j = 1, totaom(k)
            if (aom2tnk(k,j) .eq. i) then
               proceed = .true.
               if (.not. sdlfsedone(k)) then
                  calchess = .true.
               endif
               exit
            endif
         enddo
         if (proceed) then
c
c           compute Hessian elements if needed
c
            if (calchess) then
               call lfsehess (k)
               sdlfsedone(k) = .true.
            endif
c
c           increment diagonal and non-diagonal Hessian elements
c
            ia = i
            do ca = 1, 3
               do cb = 1, 3
                  do l = 1, totaom(k)
                     ib = aom2tnk(k,l)
                     call numeral (ia,pia,size)
                     call numeral (ib,pib,size)
                     call numeral (ca,pca,size)
                     call numeral (cb,pcb,size)
                     pt = pia//pib//pca//pcb
                     pt2 = pib//pia//pcb//pca
                     found = .false.
                     do m = 1, maxhlfse
                        if (hidlfse(k,m) .eq. blank16) exit
                        if ((hidlfse(k,m) .eq. pt) 
     &                      .or. (hidlfse(k,m) .eq. pt2)) then
                           if (ca .eq. 1) then
                              hessx(cb,ib) = hessx(cb,ib) + hellfse(k,m)
                           else if (ca .eq. 2) then
                              hessy(cb,ib) = hessy(cb,ib) + hellfse(k,m)
                           else if (ca .eq. 3) then
                              hessz(cb,ib) = hessz(cb,ib) + hellfse(k,m)
                           endif
                           found = .true.
                           exit
                        endif
                     enddo
                     if (.not. found) then
  10                    format (/,' LFSE4LFMM2A  --  unable to find ',
     &                          'Hessian element for LFSE ',/,
     &                          'atom1: ',i6,' coord1: ',i1,
     &                          ' atom2 : ',i6,' coord2: ',i1)
                        write(iout,10) ia,ca,ib,cb
                        call fatal
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
      return
      end subroutine
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine lfse4lfmm2b  --  slow atom-by-atom ligand field Hessian  ##   
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "lfse4lfmm2b" calculates second derivatives of the ligand
c     field stabilization energy with respect to Cartesian coordinates 
c     using finite difference method. This subroutine imply six
c     calls of ICCG routines and is more time consuming.
c     
c     
      subroutine lfse4lfmm2b (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'deriv.i'
      include 'hessn.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      integer i,j,k,l
      integer ii
      real*8 e
      real*8 stp
      real*8 old
      real*8, allocatable :: d0(:,:)
      real*8, allocatable :: d1(:,:)
      logical proceed
c
c
c     set stepsize for derivatives
c
      stp = 1.0d-5
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (d0(3,n))
      allocate (d1(3,n))

      do k = 1, nlfse
c
c        skip d0 and d10
c
         if (nelelfse(k).eq.0 .or. nelelfse(k).eq.10) then
            cycle
         endif
c
c        decide whether to compute the current term
c
         proceed = .false.
         do j = 1, totaom(k)
            if (aom2tnk(k,j) .eq. i) then
               proceed = .true.
               exit
            endif
         enddo
c
c     find numerical x-components via perturbed structure
c
         if (proceed) then
            old = x(i)
            x(i) = x(i) - 0.05d0*stp
            call lfsenggrad (k,e,d0)
            x(i) = x(i) + stp
            call lfsenggrad (k,e,d1)
            x(i) = old
            do j = 1, 3
               do l = 1, totaom(k)
                  ii = aom2tnk(k,l)
                  hessx(j,ii) = hessx(j,ii) + 
     &                          (d1(j,ii) - d0(j,ii)) / stp
               enddo
            enddo
c
c     find numerical y-components via perturbed structure
c
            old = y(i)
            y(i) = y(i) - 0.05d0*stp
            call lfsenggrad (k,e,d0)
            y(i) = y(i) + stp
            call lfsenggrad (k,e,d1)
            y(i) = old
            do j = 1, 3
               do l = 1, totaom(k)
                  ii = aom2tnk(k,l)
                  hessy(j,ii) = hessy(j,ii) + 
     &                          (d1(j,ii) - d0(j,ii)) / stp
               enddo
            enddo
c
c     find numerical z-components via perturbed structure
c
            old = z(i)
            z(i) = z(i) - 0.05d0*stp
            call lfsenggrad (k,e,d0)
            z(i) = z(i) + stp
            call lfsenggrad (k,e,d1)
            z(i) = old
            do j = 1, 3
               do l = 1, totaom(k)
                  ii = aom2tnk(k,l)
                  hessz(j,ii) = hessz(j,ii) + 
     &                          (d1(j,ii) - d0(j,ii)) / stp
               enddo
            enddo
         endif   
      enddo
c
c     perform deallocation of some local arrays
c
      deallocate (d0)
      deallocate (d1)
      return
      end subroutine
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine lfsehess  --  collect lfse contributes to the Hessian  ##   
c     ##                                                                    ##
c     ########################################################################
c
c
c     "lfsehess" collects all the contributions to the Hessian from
c     the ligand field stabilization energy term.
c
c
      subroutine lfsehess (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      include 'lfsehess.i'
      include 'iounit.i'
      integer i,j
      integer ia,ib,ca,cb
      integer ndd,idd,size
      real*8  dd
      real*8 energy
      real*8 derivs(3,n)
      real*8, allocatable :: secderivs(:)
      character*4 pia,pib,pca,pcb
      character*16 blank16
      character*16 pt
c
c
c     initialize storage of lfse Hessian
c
      blank16 = '                '
      do j = 1, maxhlfse
         hellfse(i,j) = 0.0d0
         hidlfse(i,j) = blank16
      enddo
c
c     calculated number of Hessian entries to modify
c
      ndd = 0
      do j = 1, 3*totaom(i)
         ndd = ndd + j
      enddo
c
c     allocate vector
c
      allocate(secderivs(ndd))
c
c     set flagsfor calculation of derivs
c
      doderlfse = 1
      scnderlfse = 1
c
c     use LFSE from ICCG group
c
      call runiccglfse (i,energy,derivs,secderivs)
c
c     read contributions to Hessian from single vector
c
      size = 4
      ia = 0
      ib = 0
      ca = 3
      cb = 3
      idd = 0
      do j = 1, ndd
         if (ia.eq.ib .and. ca.le.cb .and. ca.eq.3) then
            ca = 1
            ia = ia + 1 
            cb = 1
            ib = 1
         else if (ia.eq.ib .and. ca.le.cb .and. ca.lt.3) then
            ca = ca + 1 
            cb = 1
            ib = 1
         else if (cb .eq. 3) then
            cb = 1
            ib = ib + 1 
         else
            cb = cb + 1 
         end if
         call numeral (aom2tnk(i,ia),pia,size)
         call numeral (aom2tnk(i,ib),pib,size)
         call numeral (ca,pca,size)
         call numeral (cb,pcb,size)
         if (aom2tnk(i,ia) .lt. aom2tnk(i,ib)) then
            pt = pia//pib//pca//pcb
         else
            pt = pib//pia//pcb//pca
         endif
         dd = secderivs (j)

         if (dd.ne.dd) then
  40        format (/,' LFSEHESS  --  LFSE Hessian element is ',
     &                'is NaN; Check parameters for calculation',
     &                'of second derivatives in LFSE.')
            write(iout,40)
            call fatal
         endif
         idd = idd + 1
         hellfse(i,idd) = dd
         hidlfse(i,idd) = pt
      enddo 
c
c     deallocate vector
c
      deallocate(secderivs)

      return
      end subroutine
