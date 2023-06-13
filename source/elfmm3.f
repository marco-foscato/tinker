c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine elfmm3  -- calculates the energy for all LFMM terms    ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "elfmm" calculates the energy with respect to Cartesian coordinates 
c     for all additional LFMM terms and partitions the energy 
c     among the atoms.
c
c     literature references:
c
c     V. J. Burton, R. J. Deeth, C. M. Kemp and P.J. Gilbert,
c     "Molecular Mechanics for Coordination Complexes: The Impact
c     of Adding d-Electron Stabilization Energies", Journal of
c     American Chemical Society, 117, 8207-8415 (1995)
c
c
      subroutine elfmm3
      implicit none
c
c     Calculate all Metal-Ligand Stretch term
c
      call harm4lfmm3
      call morse4lfmm3
c
c     Calculate all Ligand-Ligand Repulsion term
c
      call vdw4lfmm3
      call ll4lfmm3
c
c     calculate electron pairing energy
c
      call pair4lfmm3
c
c     Calculate all Ligand Filed Stabilization terms
c
      call lfse4lfmm3
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine harm4lfmm3  --  harmonic metal-ligand bond stretching terms  ##   
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "harm4lfmm3" calculates energy for all LFMM harmonic M-L terms and partitions 
c     the energy among the atoms.
c
c
      subroutine harm4lfmm3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,ia,ib
      real*8 ideal, force
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
      logical header,huge

c
c     zero out the bond energy
c
      nehalf = 0
      ehalf = 0.0d0
      do i = 1, n
         aehalf(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the bond stretch energy
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
c
c     increment the total harmonic-LFMM energy and partition between atoms
c
            nehalf = nehalf + 1
            ehalf = ehalf + e
            aehalf(ia) = aehalf(ia) + 0.5d0*e
            aehalf(ib) = aehalf(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual LFMM Harmonic M-L Stretching',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',22x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' M-L',7x,2(i7,'-',a3),13x,2f10.4,f12.4)
            end if
         end if
      end do
      return
      end subroutine
c
c
c     #####################################################################################
c     ##                                                                                 ##
c     ##  subroutine morse4lfmm3  --  Morse metal-ligand bond stretching terms for LFMM  ##
c     ##                                                                                 ##
c     #####################################################################################
c
c
c     "morse4lfmm3" calculates the energy for all LFMM Morse M-L terms and partitions 
c     the energy among the atoms.
c
c
      subroutine morse4lfmm3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,ia,ib
      real*8 ideal,bde,alpha
      real*8 e,dt,expterm
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 lw,up,lwc,upc
      real*8 lwc3,upc3
      logical proceed
      logical header,huge
c
c     zero out the bond energy
c
      nemolf = 0
      emolf = 0.0d0
      do i = 1, n
         aemolf(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the bond stretch energy
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
c
c     add restrains
c
            lw = ideal - lowmolf
            up = ideal + upmolf
            lwc = lw*lw - rab2
            upc = rab2 - up*up
            lwc3 = lwc*lwc*lwc
            upc3 = upc*upc*upc
            if (lwc3 .gt. 0.0d0) then
              e = e + w8tmolf*lwc3
            endif
            if (upc3 .gt. 0.0d0) then
              e = e + w8tmolf*upc3
            endif
c
c     increment the total bond energy and partition over atoms
c
            nemolf = nemolf + 1
            emolf = emolf + e
            aemolf(ia) = aemolf(ia) + 0.5d0*e
            aemolf(ib) = aemolf(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual LFMM Morse M-L Stretching',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',22x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' M-L',7x,2(i7,'-',a3),13x,2f10.4,f12.4)
            end if
         end if
      end do
      return
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine vdw4lfmm3  -- van der Waals-like ligand-ligand interaction   ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "vdw4lfmm3" calculates the energy of all LFMM ligand-ligand
c     interaction terms and partitions the energy among the atoms.
c
c
      subroutine vdw4lfmm3
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,ia,ib,m
      real*8 a,b
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
      logical header,huge
c
c     zero out the energy
c
      nevwlf = 0
      evwlf = 0.0d0
      do i = 1, n
         aevwlf(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the energy
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
c
c     increment the total l-l energy and partition between atoms
c
            evwlf = evwlf + e
            nevwlf = nevwlf + 1
            aevwlf(ia) = aevwlf(ia) + 0.5d0*e
            aevwlf(ib) = aevwlf(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual LFMM vdW L-L Stretching',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',20x,
     &                       'Distance',5x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),rab,e
   20          format (' L-L',7x,2(i7,'-',a3),13x,f10.4,f12.4)
            end if
         end if
      end do
      return
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine ll4lfmm3  --  pure repulsive ligand-ligand interaction       ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "ll4lfmm3" calculates the energy of all LFMM ligand-ligand pure
c     repulsive terms and partitions the energy among the atoms.
c
c
      subroutine ll4lfmm3
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,ia,ib
      real*8 aaa,enne
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
      logical header,huge
c
c     zero out the energy
c
      nelllf = 0
      elllf = 0.0d0
      do i = 1, n
         aelllf(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the energy
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
c
c     increment the total l-l energy and partition between atoms
c
            nelllf = nelllf + 1
            elllf = elllf + e
            aelllf(ia) = aelllf(ia) + 0.5d0*e
            aelllf(ib) = aelllf(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual LFMM pure repulsive L-L ',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',20x,
     &                       'Distance',5x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),rab,e
   20          format (' L-L',7x,2(i7,'-',a3),13x,f10.4,f12.4)
            end if
         end if
      end do
      return
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine pair4lfmm3  --  empirical rule for electron pairing energy   ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "pair4lfmm3" calculate energy for all LFMM empirical electron
c     pairing terms and partitions the energy among the atoms.
c
c
      subroutine pair4lfmm3
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,ia,ib
      real*8 xab,yab,zab,dist
      real*8 e,w
      real*8 eprm(7)
      logical proceed
      logical header,huge
c
c     zero out the energy and partitioning terms
c
      nepair = 0
      epair = 0.0d0
      do i = 1, n
         aepair(i) = 0.0d0
      end do
      header = .true.
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
c     compute energy from polynomial and first derivative
c
            e = w*eprm(1) + w*eprm(2)*dist
            do j = 3, 7
               e = e + w*eprm(j)*(dist**(1-j))
            enddo
c
c     increment the total pairing energy and partition between atoms
c
            nepair = nepair + 1
            epair = epair + e
            aepair(ia) = aepair(ia) + 0.5d0*e
            aepair(ib) = aepair(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (abs(e) .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual LFMM empirical pairing ',
     &                       'energy :',
     &                    //,' Type',14x,'Atom Names',21x,
     &                       'M-L bond',4x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),dist,e
   20          format (' pair',6x,2(i7,'-',a3),13x,f10.4,f12.4)
            end if
         end if
      end do
      return
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine lfse4lfmm3  --  ligand field stabilization energy            ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "lfse4lfmm3" calculates Ligand Field Stabilization Energy for all 
c     LFMM centers and partitions the energy among the atoms 
c     (100% to central atom = the metal).
c     
c
      subroutine lfse4lfmm3
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'lfmm.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'inform.i'
      include 'energi.i'
      integer i,im
      real*8 e
      real*8 lfsenergy
      logical header,huge
c
c     zero out the energy
c
      nelfse = 0
      elfse = 0.0d0
      do i = 1, n
         aelfse(i) = 0.0d0
      enddo
      header = .true.
c
c     calculate the ligand field stabilization energy
c
      do i = 1, nlfse
         im = imlfse(i)
c
c        get ligand fiels stabilization energy
c
         e = lfsenergy (i)
c
c        increment total lfse and partition
c
         nelfse = nelfse + 1
         elfse = elfse + e
         aelfse(im) = aelfse(im) + e
c
c        print a message if the energy of this interaction is large
c
         huge = (abs(e) .gt. 5.0d0)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
                write (iout,10)
  10           format (/,' Individual LFSE terms :',
     &                //,16x,'Atom Name',16x,'Energy'/)
            endif
            write (iout,20)  im,name(im),e
  20        format (' LFSE',6x,i7,'-',a3,13x,f12.4)
         endif 
      enddo
      return
      end subroutine
