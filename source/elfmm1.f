c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine elfmm1  --  energy & derivatives for all LFMM terms  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "elfmm1" calls all the routines calculating the energy and 
c     the first derivatives with respect to Cartesian coordinates 
c     for all LFMM terms.
c
c     literature references:
c
c     V. J. Burton, R. J. Deeth, C. M. Kemp and P.J. Gilbert,
c     "Molecular Mechanics for Coordination Complexes: The Impact
c     of Adding d-Electron Stabilization Energies", Journal of
c     American Chemical Society, 117, 8207-8415 (1995)
c
c
      subroutine elfmm1
      implicit none
c
c     Calculate all Metal-Ligand Stretch term
c
      call harm4lfmm1
      call morse4lfmm1
c
c     Calculate all Ligand-Ligand Repulsion term
c
      call vdw4lfmm1
      call ll4lfmm1
c
c     calculate electron pairing energy
c
      call pair4lfmm1
c
c     Calculate all Ligand Filed Stabilization terms
c
      call lfse4lfmm1
      end subroutine
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine harm4lfmm1  --  LFMM harmonic M-L stretch & derivatives  ##   
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "harm4lfmm1" calculates energy and first derivatives with respect
c     to Cartesian coordinates for all LFMM harmonic M-L terms, and 
c     increments the energy, first derivatives and virial.
c
c
      subroutine harm4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib
      real*8 ideal, force
      real*8 e,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the bond energy
c
      ehalf = 0.0d0
      do i = 1, n
         dehalf(1,i) = 0.0d0
         dehalf(2,i) = 0.0d0
         dehalf(3,i) = 0.0d0
      end do
c
c     calculate the bond stretch energy and first derivatives
c     
      do i = 1, nhalf
         ia = ihalf(1,i)
         ib = ihalf(2,i)
         ideal = lhalf(i)
         force = khalf(i)
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
            rab = sqrt(xab*xab + yab*yab + zab*zab)
c
c     Harmonic interaction has functional form U(r) = 1/2 * k * (r - r0)**2
c
            e = 0.5d0*force*(rab-ideal)**2
            deddt = force*(rab-ideal)
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total harmonic-LFMM energy and first derivatives
c
            ehalf = ehalf + e
            dehalf(1,ia) = dehalf(1,ia) + dedx
            dehalf(2,ia) = dehalf(2,ia) + dedy
            dehalf(3,ia) = dehalf(3,ia) + dedz
            dehalf(1,ib) = dehalf(1,ib) - dedx
            dehalf(2,ib) = dehalf(2,ib) - dedy
            dehalf(3,ib) = dehalf(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end subroutine
c
c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine morse4lfmm1  --  LFMM Morse M-L stretch energy  & derivatives  ##
c     ##                                                                            ##
c     ################################################################################
c
c
c     "morse4lfmm1" calculates the energy and first derivatives with respect 
c     to Cartesian coordinates for all LFMM Morse Metal-Ligand bond stretch 
c     terms and increments the energy, first derivatives and virial.
c
c
      subroutine morse4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'lfmm.i'
      integer i,ia,ib
      real*8 ideal,bde,alpha
      real*8 e,expterm
      real*8 dt,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 rab,rab2
      real*8 lw,up,lwc,upc
      real*8 lwc2,lwc3,upc2,upc3
      logical proceed
c
c
c     zero out the bond energy
c
      emolf = 0.0d0
      do i = 1, n
         demolf(1,i) = 0.0d0
         demolf(2,i) = 0.0d0
         demolf(3,i) = 0.0d0
      end do
c
c     calculate the bond stretch energy and first derivatives
c
      do i = 1, nmolf
         ia = imolf(1,i)
         ib = imolf(2,i)
         ideal = lmolf(i)
         bde = bmolf(i)
         alpha = amolf(i)
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
            e = bde*((1.0d0-expterm)**2-1.0d0) 
            deddt = 2.0d0*bde*alpha*(1.0d0-expterm)*expterm
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
              e = e + w8tmolf*lwc3
              deddt = deddt - 6.0d0*w8tmolf*rab*lwc2
            endif
            if (upc3 .gt. 0.0d0) then
              e = e + w8tmolf*upc3
              deddt = deddt + 6.0d0*w8tmolf*rab*upc2
            endif
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total bond energy and first derivatives
c
            emolf = emolf + e
            demolf(1,ia) = demolf(1,ia) + dedx
            demolf(2,ia) = demolf(2,ia) + dedy
            demolf(3,ia) = demolf(3,ia) + dedz
            demolf(1,ib) = demolf(1,ib) - dedx
            demolf(2,ib) = demolf(2,ib) - dedy
            demolf(3,ib) = demolf(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end subroutine
c
c
c     ###############################################################################
c     ##                                                                           ##
c     ##  subroutine vdw4lfmm1  --  LFMM vdW-like L-L interaction and derivatives  ##
c     ##                                                                           ##
c     ###############################################################################
c
c
c     "vdw4lfmm1" calculates the energy and first derivatives with respect
c     to Cartesian coordinates for all LFMM ligand-ligand interactions term
c     and increments the energy, first derivatives and virial.
c
c
      subroutine vdw4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'lfmm.i'
      integer i,ia,ib,m
      real*8 a,b
      real*8 e,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the energy
c
      evwlf = 0.0d0
      do i = 1, n
         devwlf(1,i) = 0.0d0
         devwlf(2,i) = 0.0d0
         devwlf(3,i) = 0.0d0
      end do
c
c     calculate the energy and first derivatives
c  
      do i = 1, nvwlf
         ia = ivwlf(1,i)
         ib = ivwlf(2,i)
         m = mvwlf(i)
         a = avwlf(i)
         b = bvwlf(i)
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
            rab = sqrt(xab*xab + yab*yab + zab*zab)
c
c     Ligand-Ligand interaction as Lennard-Jones-like term
c     U(r) = a/r**m-b/r**6
c
            e = (a/(rab**m))-(b/(rab**6))
            deddt = -a*m*rab**(-m-1) + 6*b*(1/rab**7)
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total l-l energy and first derivatives
c
            evwlf = evwlf + e
            devwlf(1,ia) = devwlf(1,ia) + dedx
            devwlf(2,ia) = devwlf(2,ia) + dedy
            devwlf(3,ia) = devwlf(3,ia) + dedz
            devwlf(1,ib) = devwlf(1,ib) - dedx
            devwlf(2,ib) = devwlf(2,ib) - dedy
            devwlf(3,ib) = devwlf(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end subroutine
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine ll4lfmm1  --  LFMM pure L-L repulsion & derivatives  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "ll4lfmm1" calculates energy and first derivatives with respect
c     to Cartesian coordinates for all LFMM ligand-ligand pure repulsive
c     terms and increments the energy, first derivatives and virial.
c
c
      subroutine ll4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'lfmm.i'
      integer i,ia,ib
      real*8 aaa,enne
      real*8 e,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the energy
c
      elllf = 0.0d0
      do i = 1, n
         delllf(1,i) = 0.0d0
         delllf(2,i) = 0.0d0
         delllf(3,i) = 0.0d0
      end do
c
c     calculate the energy and first derivatives
c     
      do i = 1, nlllf
         ia = illlf(1,i)
         ib = illlf(2,i)
         enne = mlllf(i)
         aaa = alllf(i)
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
            rab = sqrt(xab*xab + yab*yab + zab*zab)
c
c     Ligand-Ligand interaction has functional form U(r) = A*r**(-n) .
c
            e = aaa*rab**(-enne)
            deddt = -enne*aaa*rab**(-enne-1)
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total l-l energy and first derivatives
c
            elllf = elllf + e
            delllf(1,ia) = delllf(1,ia) + dedx
            delllf(2,ia) = delllf(2,ia) + dedy
            delllf(3,ia) = delllf(3,ia) + dedz
            delllf(1,ib) = delllf(1,ib) - dedx
            delllf(2,ib) = delllf(2,ib) - dedy
            delllf(3,ib) = delllf(3,ib) - dedz
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end subroutine
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine pair4lfmm1  --  LFMM  electron pairing and derivatives  ##
c     ##                                                                     ##
c     #########################################################################
c
c
c     "pair4lfmm1" calculates the energy and first derivatives with respect
c     to Cartesian coordinates for all LFMM empirical electron pairing 
c     terms and increments the energy, first derivatives and virial.
c
c
      subroutine pair4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'lfmm.i'
      integer i,j,ia,ib
      real*8 xab,yab,zab,dist
      real*8 e,w,deddist,de
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 eprm(7)
      logical proceed
c
c
c     zero out the energy and first derivatives
c
      epair = 0.0d0
      do i = 1, n
         depair(1,i) = 0.0d0
         depair(2,i) = 0.0d0
         depair(3,i) = 0.0d0
      end do
c
c     calculate the electron pairing energy
c
      do i = 1, npair
         ia = ipair(1,i)
         ib = ipair(2,i)
         eprm = ppair(i,:)
         w = wpair(i)
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
            dist = sqrt(xab*xab + yab*yab + zab*zab)
c
c     compute energy from polynomial and first derivatives
c
            e = w*eprm(1) + w*eprm(2)*dist
            deddist = w*eprm(2)
            do j = 3, 7
               e = e + w*eprm(j)*(dist**(1-j))
               deddist = deddist - w*((j-1)*eprm(j)*(dist**(-j)))
            enddo
c
c     compute chain rule terms needed for derivatives
c
            if (dist .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddist / dist
            endif
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total pairing energy and first derivatives
c
            epair = epair + e
            depair(1,ia) = depair(1,ia) + dedx
            depair(2,ia) = depair(2,ia) + dedy
            depair(3,ia) = depair(3,ia) + dedz
            depair(1,ib) = depair(1,ib) - dedx
            depair(2,ib) = depair(2,ib) - dedy
            depair(3,ib) = depair(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end subroutine
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine lfse4lfmm1  --  LFSE and derivatives  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "lfse4lfmm1" calculates energy and first derivatives of the
c     Ligand Field Stabilization Energy with
c     respect to Cartesian coordinates for all LFMM centers,
c     and increments the energy, first derivatives and virial.
c
c
      subroutine lfse4lfmm1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'deriv.i'
      include 'virial.i'
      include 'lfmm.i'
      integer i,j,k
      real*8 e
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 g(3,n)
c
c
c     zero out the energy and derivatives
c
      elfse = 0.0d0
      do i = 1, n
         do j = 1, 3
            delfse(j,i) = 0.0d0
            g(j,i) = 0.0d0
         enddo
      enddo
c
c     calculate the lfse and derivatives for all LFMM centers
c
      do i = 1, nlfse
         call lfsenggrad (i,e,g)
c
c        increment the total lfse energy and derivatives
c
         elfse = elfse + e
         do j = 1, totaom(i)
            do k = 1, 3
               delfse(k,aom2tnk(i,j)) = delfse(k,aom2tnk(i,j)) 
     &                                  + g(k,aom2tnk(i,j))
            enddo
         enddo
c
c        increment the internal virial tensor components
c
         do j = 1, totaom(i)
            vxx = x(aom2tnk(i,j)) * g(1,aom2tnk(i,j))
            vyx = y(aom2tnk(i,j)) * g(1,aom2tnk(i,j))
            vzx = z(aom2tnk(i,j)) * g(1,aom2tnk(i,j))
            vyy = y(aom2tnk(i,j)) * g(2,aom2tnk(i,j))
            vzy = z(aom2tnk(i,j)) * g(2,aom2tnk(i,j))
            vzz = z(aom2tnk(i,j)) * g(3,aom2tnk(i,j))
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         enddo
      enddo
      return
      end subroutine
c
c
c     ###################################################################################
c     ##                                                                               ##
c     ##  subroutine lfsenggrad  --  energy and gradient of ligand field stab. energy  ##
c     ##                                                                               ##
c     ###################################################################################
c
c
c     "lfsenggrad" collects the ligand field stabilization energy (according
c     to LFMM method) and first derivatives with respect to Cartesian
c     coordinates for a single LFMM center.
c
c
      subroutine lfsenggrad (i,e,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      integer i,j,k
      real*8 e
      real*8 g(3,n)
      real*8 energy
      real*8, allocatable, dimension(:,:) :: derivs
      real*8, allocatable, dimension(:) :: secderivs
c
c
c     zero out the energy and grad
c
      e = 0.0d0
      do k = 1, n
         do j = 1, 3
            g(j,k) = 0.0d0
         end do
      end do
c
c     skip d0 and d10
c
      if (nelelfse(i).eq.0 .or. nelelfse(i).eq.10) then
         return
      endif
c
c     allocate used and dummy arrays
c
      allocate(derivs(3,totaom(i)))
      allocate(secderivs(1))
c
c     must set flag for calculation of derivs
c
      doderlfse = 1
      scnderlfse = 0
c
c     zero out the energy and grad
c
      do k = 1, totaom(i)
         do j = 1, 3
            derivs(j,k) = 0.0d0
         end do
      end do
c
c     use LFSE from ICCG group
c
      call runiccglfse (i,energy,derivs,secderivs)
      e = energy
c
c     reorder derivatives
c
      do j = 1,totaom(i)
        g(1,aom2tnk(i,j)) = derivs(1,j)
        g(2,aom2tnk(i,j)) = derivs(2,j)
        g(3,aom2tnk(i,j)) = derivs(3,j)
      enddo
c
c     deallocate
c
      deallocate(derivs)
      deallocate(secderivs)

      return
      end subroutine
