c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
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
c     "iccglfse" is a routinge acting as a bridge between Tinker code
c     and LFSE that calculates the ligand field stabilization energy
c     and derivatives and is developed at ICCG by Rob J. Deeth and
c     collaborators. 
c     This file is only an empty cover. You can obtain the source code
c     from Rob J. Deeth (University of Warwick) at r.j.deeth@warwick.ac.uk
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
     &                     dbiccg,
     &                     stpfindif,
     &                     dispfindif,
     &                     secsderivat,
     &                     stabeng,
     &                     numdelectrons,
     &                     spinmultipl,
     &                     evaleigval,
     &                     energy,
     &                     derivs,
     &                     hessvec,
     &                     aspre2014)
      implicit none
      include "iounit.i"
      integer totaom,nligands,numbari
      integer dispfindif,secsderivat,stabeng
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
      logical sjflag,adflag,dbiccg,doderivflag,aspre2014

      write(iout,*) ' '
      write(iout,*) ' ICCGLFSE -- Please obtain ICCGLFSE source code'
      write(iout,*) ' from Rob J. Deeth (University of Warwick)'
      write(iout,*) ' and recompile Tinker before attempting LFMM'
      write(iout,*) ' calculations.'
      write(iout,*) ' '
      write(iout,*) ' Write to r.j.deeth@warwick.ac.uk'
      write(iout,*) ' '

      call fatal

      end subroutine
