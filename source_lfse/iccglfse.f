c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Rob J. Deeth        ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine iccglfse  --  ligand field stabilization energy by ICCG  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "iccglfse" is a routine acting as a bridge between Tinker code
c     the code developed at ICCG by Rob J. Deeth and co.
c     for the calculation of the ligand field stabilization energy LFSE
c     and its derivatives.
c     This code define the subroutine 'iccglfse' that overwrites the 
c     corresponding, yet empty (not calculating anything) 
c     subroutine defined in the
c     code distributed within the Tinker package.
c
c
      subroutine iccglfse (totaom,
     &                     coordx, 
     &                     coordy,
     &                     coordz,
     &                     nligands,
     &                     idlig,
     &                     typlig,
     &                     neiglig,
     &                     nillig,
     &                     nsslig,
     &                     lpid2,
     &                     lpid3,
     &                     chrglig,
     &                     aomsig,
     &                     aompix,
     &                     aompiy,
     &                     aomxds,
     &                     doderivflag,
     &                     sjflag,
     &                     adflag,
     &                     numbari,
     &                     dbugflag,
     &                     stpfindif,
     &                     dispfindif,
     &                     secderivat,
     &                     stabeng,
     &                     numdelectrons,
     &                     spinmultipl,
     &                     evaleigval,
     &                     energy,
     &                     derivs,
     &                     hessvec,
     &                     aspre2014)
      implicit double precision (a-h,o-z)
      include 'units.ins'
      include 'param.ins'
      include 'switch.ins'
      include 'numinf.ins'
      include 'clf.ins'
      include 'secd.ins'
      include 'eigv.ins'
      include 'anal1d.ins'
      integer ii
      integer totaom,nligands,numbari
      integer dispfindif,secderivat,stabeng
      integer numdelectrons,spinmultipl,evaleigval
      integer idlig(*),typlig(*)
      integer neiglig(*),nillig(*),nsslig(*)
      integer lpid2(*),lpid3(*)
      real*8 energy
      real*8 stpfindif
      real*8 coordx(*),coordy(*),coordz(*)
      real*8 chrglig(*)
      real*8 hessvec(*)
      real*8 derivs(3,*)
      real*8 aomsig(7,*),aompix(7,*),aompiy(7,*),aomxds(7,*)
      real*8 xyzout(3,totaom)
      logical sjflag,adflag,dbugflag,doderivflag,aspre2014
c
c
c     reproduce old behaviour
c
      reproduceold=aspre2014
c
c     number of atoms in AOM region
c
      natom=totaom
c
c     coordinates of re-numbered atom list 
c
      do ii=1, natom
         xyzout(1,ii) = coordx(ii)
         xyzout(2,ii) = coordy(ii)
         xyzout(3,ii) = coordz(ii)
      enddo
c
c     number of coordinating atoms
c
      nligs=nligands
c
c     properties of ligands
c
      do ii=1, nligs
        ligdef(ii,1) = typlig(ii)
        ligdef(ii,2) = idlig(ii)
        ligdef(ii,3) = neiglig(ii)
        ligdef(ii,4) = nillig(ii)
        ligdef(ii,5) = nsslig(ii)
      enddo
c
c     count ligand types
c
      nltyp = 0
      do ii=1, nligs
        if (typlig(ii) .gt. nltyp) nltyp = typlig(ii)
      enddo
c
c     subsidiary atoms for ligands
c
      nslugs = 0
      slgdef = 0 
      nslugp = 0
      do n=1, nligs
         if ((nsslig(n).eq.1) .or. (nsslig(n) .eq. 2)) then
            nslugp = nslugp - 1
            nslugs             = nslugs - nsslig(n)
            slgdef(n,1)        = lpid2(n)
            slgdef(n,3)        = nslugp
            slug(nslugp,1)     = lpid2(n)
            slug(nslugp,2)     = n
            if (nsslig(n) .eq. 2) then
               nslugp = nslugp - 1
               slgdef(n,2)        = lpid3(n)
               slgdef(n,4)        = nslugp
               slug(nslugp,1)     = lpid3(n)
               slug(nslugp,2)     = n
            endif
         endif
      enddo
c
c     clf charges on metal and ligands
c
      do ii=1, nligs+1
         clfchg(ii)=chrglig(ii)
      enddo
c
c     aom parameters
c      
      do ii=1, 7
         do ij=1, nltyp 
            esig(ij,ii) = aomsig(ii,ij)
            epix(ij,ii) = aompix(ii,ij)
            epiy(ij,ii) = aompiy(ii,ij)
            exds(ij,ii) = aomxds(ii,ij)
         enddo
      enddo
c
c     Calculate derivatives flag
c
      doderivs=doderivflag
c
c     Shaffer and Jorgensen d/s-mixing flag
c
      schjor=sjflag
c
c     analytic derivative flag
c
      an1std=adflag
c
c     baricenter of Shaffer and Jorgensen d/s-mixing
c
      nbarry=numbari
c
c     debug flag
c
      dbug=dbugflag
c
c     step size for finite difference
c
      drvstp=stpfindif
c
c     displacement for finite difference
c
      nderiv=dispfindif
c
c     controls derivation steps (1st and 2nd derivs)
c
      nwiccg=secderivat
c
c     controls stabilizzation energy formula
c
      nclfse=stabeng
c
c     number of d-elecrons
c
      ndelec=numdelectrons
c
c     spin multiplicity
c
      mult=spinmultipl
c
c     evaluation eigenvalues
c
      evaleig=evaleigval
c
c     initialize results
c
      energy = 0.0d0
      dd = 0.0d0
c
c     finally submit calculation to CLFSE routine
c

      call clfse(xyzout,natom,derivs,energy)

c
c     collect second derivatives
c
      if (nwiccg .eq. 1) then
         ndd = 0
         do ii=1, 3*natom
            ndd = ndd + ii
         enddo
         do ii=1, ndd
            hessvec(ii) = dd(ii)
         enddo
      endif

      return
      end subroutine
