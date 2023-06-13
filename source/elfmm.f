c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine elfmm  -- calculates the energy for all LFMM terms  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "elfmm" calls all the routines that calculate the energy for all 
c     LFMM terms.
c
c     Literature references:
c
c     V. J. Burton, R. J. Deeth, C. M. Kemp and P.J. Gilbert,
c     "Molecular Mechanics for Coordination Complexes: The Impact
c     of Adding d-Electron Stabilization Energies", Journal of
c     American Chemical Society, 117, 8207-8415 (1995)
c
c
      subroutine elfmm
      implicit none
c
c     Calculate all Metal-Ligand Stretch term
c
      call harm4lfmm
      call morse4lfmm
c
c     Calculate all Ligand-Ligand Repulsion term
c
      call vdw4lfmm
      call ll4lfmm
c
c     calculate electron pairing energy
c
      call pair4lfmm
c
c     Calculate all Ligand Filed Stabilization terms
c
      call lfse4lfmm
      end subroutine
c
c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine harm4lfmm  --  harmonic metal-ligand stretch potential energy  ##   
c     ##                                                                            ##
c     ################################################################################
c
c
c     "harm4lfmm" calculates energy for all LFMM harmonic metal-ligand bond
c     stretching terms.
c
c
      subroutine harm4lfmm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      integer i,ia,ib
      real*8 ideal, force
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the bond energy
c
      ehalf = 0.0d0
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
            ehalf = ehalf + e
         end if
      end do
      return
      end subroutine
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine morse4lfmm  --  Morse metal-ligand stretch potential energy  ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "morse4lfmm" calculates the energy for all LFMM Morse metal-ligand bond
c     stretching terms.
c
c
      subroutine morse4lfmm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      integer i,ia,ib
      real*8 ideal,bde,alpha
      real*8 e,dt,expterm
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 lw,up,lwc,upc
      real*8 lwc3,upc3
      logical proceed
c
c
c     zero out the bond energy
c
      emolf = 0.0d0
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
c     increment the total bond energy
c
            emolf = emolf + e
         end if
      end do
      return
      end subroutine
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine vdw4lfmm  --  vdW-like ligand-ligand potential energy  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "vdw4lfmm" calculates the energy for all LFMM van der Waals-like
c     ligand-ligand interactions.
c
c
      subroutine vdw4lfmm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      integer i,ia,ib,m
      real*8 a,b
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
c
c     zero out the energy
c
      evwlf = 0.0d0
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
c     increment the total l-l energy
c
            evwlf = evwlf + e
         end if
      end do
      return
      end subroutine
c
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine ll4lfmm  --  pure repulsive ligand-ligand potential energy  ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "ll4lfmm" calculates the potential energy for LFMM lignad-ligands
c     pure repulsive terms.
c
c
      subroutine ll4lfmm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      integer i,ia,ib
      real*8 aaa,enne
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the energy
c
      elllf = 0.0d0
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
c     increment the total l-l energy
c
            elllf = elllf + e
         end if
      end do
      return
      end subroutine
c
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine pair4lfmm  --  empirical rule for electron pairing energy   ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "pair4lfmm" calculate potential energy for all LFMM empiric terms
c     accounting for the electron pairing energy.
c
c
      subroutine pair4lfmm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      include 'energi.i'
      include 'lfmm.i'
      integer i,j,ia,ib
      real*8 xab,yab,zab,dist
      real*8 e,w
      real*8 eprm(7)
      logical proceed
c
c
c     zero out the energy
c
      epair = 0.0d0
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
c     increment the total pairing energy
c
            epair = epair + e
         end if
      end do
      return
      end subroutine
c
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine lfse4lfmm  --  ligand field stabilization energy            ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "lfse4lfmm" calculates potential energy for all LFMM centers
c     
c
      subroutine lfse4lfmm
      implicit none
      include 'sizes.i'
      include 'lfmm.i'
      include 'energi.i'
      integer i
      real*8 e
      real*8 lfsenergy
c
c     zero out the energy
c
      elfse = 0.0d0
c
c     collect the ligand field stabilization energies
c
      do i = 1, nlfse
         e = lfsenergy (i)
         elfse = elfse + e
      enddo
      return
      end subroutine
c
c
c     ####################################################################
c     ##                                                                ##
c     ##   function lfsenergy  --  evaluates ligand field stab. energy  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "lfsenergy" returns the ligand field stabilization energy (according
c     to LFMM method) of a given LFMM center.
c
c
      function lfsenergy (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      integer i
      real*8 e
      real*8 lfsenergy
      real*8, allocatable, dimension(:,:) :: g
      real*8, allocatable, dimension(:) :: h
c
c
c     allocate dummy array
c
      allocate(g(3,totaom(i)))
      allocate(h(1))
c
c     initialize
c
      e = 0.0d0
      g = 0.0d0
      h = 0.0d0
c
c     must set flag for calculation of derivs to false 
c
      doderlfse = 0
      scnderlfse = 0
c
c     use LFSE from ICCG group
c
      call runiccglfse (i,e,g,h)

      lfsenergy = e
c
c     deallocate
c
      deallocate(g)
      deallocate(h)

      return
      end function
c
c
c     ###################################################################
c     ##                                                               ##
c     ##   subroutine  runiccglfse  --  run LFSE routines from ICCG    ##
c     ##                                                               ##
c     ###################################################################
c
c     "runiccglfse" works as interface between Tinker and the subroutines
c     developed by Rob J. Deeth (ICCG - University of Warwick).
c     This routine converts all the information needed to calculate the
c     lignd field stabilization energy (and derivatives) in a form 
c     compatible to the vector-based data structure in ICCG routines, 
c     then calls such routines and collects the results.
c
c
      subroutine runiccglfse (i,iccge,iccgg,iccgh)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      integer i
      integer ii,ij
      integer tnk2aom
      integer unqtype
      integer, allocatable, dimension(:) :: idlig,typlig,neiglig,nillig
      integer, allocatable, dimension(:) :: nsslig,lpid2,lpid3
      real*8 iccge
      real*8 iccgh(*)
      real*8 iccgg(3,*)
      real*8, allocatable, dimension(:) :: coordx,coordy,coordz,chrglig
      real*8, allocatable, dimension(:,:) :: aomsig,aompix,aompiy,aomxds
      logical sjflag,adflag,dbiccg,doderivflag
c
c
c     skip d0 and d10
c
      if (nelelfse(i).eq.0 .or. nelelfse(i).eq.10) then
         return
      endif
c
c     prepare vectors for LFSE routine from ICCG
c
      allocate (coordx(totaom(i)))
      allocate (coordy(totaom(i)))
      allocate (coordz(totaom(i)))

      do ii=1,3
        do ij=1,totaom(i)
          if (ii.eq.1) coordx(ij) = x(aom2tnk(i,ij))
          if (ii.eq.2) coordy(ij) = y(aom2tnk(i,ij))
          if (ii.eq.3) coordz(ij) = z(aom2tnk(i,ij))
        enddo
      enddo

      allocate (idlig(nligs(i)))
      allocate (typlig(nligs(i)))
      allocate (neiglig(nligs(i)))
      allocate (nillig(nligs(i)))
      allocate (nsslig(nligs(i)))
      allocate (lpid2(nligs(i)))
      allocate (lpid3(nligs(i)))
      allocate (chrglig(nligs(i) + 1))

      do ii = 1, nligs(i)
        idlig(ii) = tnk2aom(i,illfse(i,ii))
      enddo
      do ii = 1, nligs(i)
        typlig(ii) = ligtyp(i,ii)
      enddo
      do ii = 1, nligs(i)
        neiglig(ii) = tnk2aom(i,imlfse(i))
      enddo
      call getligsnotinline (i,nillig)
      do ii = 1, nligs(i)
        nsslig(ii) = nss1(i,ii) + nss2(i,ii)
      enddo
      do ii = 1, nligs(i)
        if (nss1(i,ii) .eq. 0) then
           lpid2(ii) = 0
        else
           lpid2(ii) = tnk2aom(i,iss1(i,ii))
        endif
      enddo
      do ii = 1, nligs(i)
         if (nss2(i,ii) .eq. 0) then
            lpid3(ii) = 0
         else
            lpid3(ii) = tnk2aom(i,iss2(i,ii))
         endif
      enddo
c
c     No need to send formal charge of metal and ligand atoms
c     because electrostatics is calculated elsewhere (see echarge)
c
      do ii = 1, nligs(i) + 1
        chrglig(ii) = 0.0d0
      enddo
      allocate (aomsig(7,nligs(i)))
      allocate (aompix(7,nligs(i)))
      allocate (aompiy(7,nligs(i)))
      allocate (aomxds(7,nligs(i)))
      do ii = 1,7
         do ij = 1, maxlig
            unqtype = ligunqtyps(i,ij,1)
            if (unqtype .lt. 0) then
               exit
            else
               aomsig(ii,ij) = esig(i,ligunqtyps(nlfse,ij,2),ii)
            endif
         enddo
      enddo
      do ii = 1,7
         do ij = 1, maxlig
            unqtype = ligunqtyps(i,ij,1)
            if (unqtype .lt. 0) then
               exit
            else
               aompix(ii,ij) = epix(i,ligunqtyps(nlfse,ij,2),ii)
            endif
         enddo
      enddo
      do ii = 1,7
         do ij = 1, maxlig
            unqtype = ligunqtyps(i,ij,1)
            if (unqtype .lt. 0) then
               exit
            else
               aompiy(ii,ij) = epiy(i,ligunqtyps(nlfse,ij,2),ii)
            endif
         enddo
      enddo
      do ii = 1,7
         do ij = 1, maxlig
            unqtype = ligunqtyps(i,ij,1)
            if (unqtype .lt. 0) then
               exit
            else
               aomxds(ii,ij) = exds(i,ligunqtyps(nlfse,ij,2),ii)
            endif
         enddo
      enddo
      doderivflag = .false.
      sjflag = .false.
      adflag = .false.
      dbiccg = .false.
      if (doderlfse .eq. 1) doderivflag = .true.
      if (shajorlfse .eq. 1 ) sjflag = .true.
      if (anlderlfse .eq. 1 ) adflag = .true.
      if (debuglfse .eq. 1 ) dbiccg = .true.
c
c     use LFSE from ICCG group
c
      call iccglfse (totaom(i),
     &               coordx,
     &               coordy,
     &               coordz,
     &               nligs(i),
     &               idlig,
     &               typlig,
     &               neiglig,
     &               nillig,
     &               nsslig,
     &               lpid2,
     &               lpid3,
     &               chrglig,
     &               aomsig,
     &               aompix,
     &               aompiy,
     &               aomxds,
     &               doderivflag,
     &               sjflag,
     &               adflag,
     &               barysjlfse,
     &               dbiccg,
     &               stpzlfse,
     &               displfse,
     &               scnderlfse,
     &               stbenglfse,
     &               nelelfse(i),
     &               spinlfmm(i),
     &               0,
     &               iccge,
     &               iccgg,
     &               iccgh,
     &               reproduceold)
c
c     deallocate local vectors
c
      deallocate (coordx)
      deallocate (coordy)
      deallocate (coordz)
      deallocate (idlig)
      deallocate (typlig)
      deallocate (neiglig)
      deallocate (nillig)
      deallocate (nsslig)
      deallocate (lpid2)
      deallocate (lpid3)
      deallocate (chrglig)
      deallocate (aomsig)
      deallocate (aompix)
      deallocate (aompiy)
      deallocate (aomxds)

      return
      end subroutine
c
c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine getligsnotinline  -- identify ligands forming non-flat angles  ##
c     ##                                                                            ##
c     ################################################################################
c
c     "getligsnotinline" prepare a vector with the index of atoms (A) 
c     such that the angles around the LFMM metal center (B), and 
c     with respect to a given list of atoms (C), are 0 < ABC < 175
c
c     idmt   index of LFMM center: this defines the central metal (B) 
c            and the list of ligands (C)
c     nfligs output list if indexes of atoms (A)
c
c
      subroutine getligsnotinline (idmt,lig2index)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'math.i'
      include 'couple.i'
      include 'lfmm.i'
      integer a,b,c
      integer i,j
      integer tnk2aom
      integer idmetal,idmt
      integer lig2index(12)
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xp,yp,zp,rp
      real*8 rab2,rcb2
      real*8 cosine,dot,angle
c
c
c     loop through Li-M-Lj combinations
c
      idmetal = imlfse(idmt)
      b = idmetal
      do i=1,n12(idmetal)
        a = i12(i,idmetal)
        do j=1,n12(idmetal)
          c = i12(j,idmetal)
c
c     calculate La-Mb-Lc angles
c
          xia = x(a)
          yia = y(a)
          zia = z(a)
          xib = x(b)
          yib = y(b)
          zib = z(b)
          xic = x(c)
          yic = y(c)
          zic = z(c)
          xab = xia - xib
          yab = yia - yib
          zab = zia - zib
          xcb = xic - xib
          ycb = yic - yib
          zcb = zic - zib
          rab2 = xab*xab + yab*yab + zab*zab
          rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
          if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            rp = max(rp,0.000001d0)
            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / sqrt(rab2*rcb2)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
c
c      get only the first non-linear combination
c
            if ((angle .gt. 0.0d0) .and. (angle .lt. 175.0d0)) then
              lig2index(i) = tnk2aom(idmt,c)
              exit
            endif
          endif
        enddo
      enddo
      end subroutine
