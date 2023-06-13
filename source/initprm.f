c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine initprm  --  initialize force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "initprm" completely initializes a force field by setting all
c     parameters to zero and using defaults for control values
c
c
      subroutine initprm
      implicit none
      include 'sizes.i'
      include 'angpot.i'
      include 'bndpot.i'
      include 'chgpot.i'
      include 'fields.i'
      include 'kanang.i'
      include 'kangs.i'
      include 'katoms.i'
      include 'kbonds.i'
      include 'kchrge.i'
      include 'kdipol.i'
      include 'khbond.i'
      include 'kiprop.i'
      include 'kitors.i'
      include 'kmulti.i'
      include 'kopbnd.i'
      include 'kopdst.i'
      include 'korbs.i'
      include 'kpitor.i'
      include 'kpolr.i'
      include 'kstbnd.i'
      include 'ksttor.i'
      include 'ktorsn.i'
      include 'ktrtor.i'
      include 'kurybr.i'
      include 'kvdws.i'
      include 'kvdwpr.i'
      include 'math.i'
      include 'merck.i'
      include 'mplpot.i'
      include 'polpot.i'
      include 'rxnpot.i'
      include 'solute.i'
      include 'urypot.i'
      include 'torpot.i'
      include 'units.i'
      include 'vdwpot.i'
      include 'klfmms.i'
      include 'ksnbnd.i'
      include 'lfmmset.i'
      include 'kgeoms.i'
      integer i,j
      character*3 blank3
      character*8 blank8
      character*12 blank12
      character*16 blank16
      character*20 blank20
      character*24 blank24
c
c
c     define blank character strings of various lengths
c
      blank3 = '   '
      blank8 = '        '
      blank12 = '            '
      blank16 = '                '
      blank20 = '                    '
      blank24 = '                        '
c
c     initialize strings of parameter atom types and classes
c
      do i = 1, maxnvp
         kvpr(i) = blank8
      end do
      do i = 1, maxnhb
         khb(i) = blank8
      end do
      do i = 1, maxnb
         kb(i) = blank8
      end do
      do i = 1, maxnb5
         kb5(i) = blank8
      end do
      do i = 1, maxnb4
         kb4(i) = blank8
      end do
      do i = 1, maxnb3
         kb3(i) = blank8
      end do
      do i = 1, maxnel
         kel(i) = blank12
      end do
      do i = 1, maxna
         ka(i) = blank12
      end do
      do i = 1, maxna5
         ka5(i) = blank12
      end do
      do i = 1, maxna4
         ka4(i) = blank12
      end do
      do i = 1, maxna3
         ka3(i) = blank12
      end do
      do i = 1, maxnaf
         kaf(i) = blank12
      end do
      do i = 1, maxnsb
         ksb(i) = blank12
      end do
      do i = 1, maxnu
         ku(i) = blank12
      end do
      do i = 1, maxnopb
         kopb(i) = blank8
      end do
      do i = 1, maxnopd
         kopd(i) = blank16
      end do
      do i = 1, maxndi
         kdi(i) = blank16
      end do
      do i = 1, maxnti
         kti(i) = blank16
      end do
      do i = 1, maxnt
         kt(i) = blank16
      end do
      do i = 1, maxnt5
         kt5(i) = blank16
      end do
      do i = 1, maxnt4
         kt4(i) = blank16
      end do
      do i = 1, maxnpt
         kpt(i) = blank8
      end do
      do i = 1, maxnbt
         kbt(i) = blank16
      end do
      do i = 1, maxntt
         ktt(i) = blank20
      end do
      do i = 1, maxnd
         kd(i) = blank8
      end do
      do i = 1, maxnd5
         kd5(i) = blank8
      end do
      do i = 1, maxnd4
         kd4(i) = blank8
      end do
      do i = 1, maxnd3
         kd3(i) = blank8
      end do
      do i = 1, maxnmp
         kmp(i) = blank12
      end do
      do i = 1, maxnpi
         kpi(i) = blank8
      end do
      do i = 1, maxnpi5
         kpi5(i) = blank8
      end do
      do i = 1, maxnpi4
         kpi4(i) = blank8
      end do
c
c     initialize values of some force field parameters
c
      forcefield = blank20
      do i = 1, maxtyp
         symbol(i) = blank3
         atmcls(i) = 0
         atmnum(i) = 0
         weight(i) = 0.0d0
         ligand(i) = 0
         describe(i) = blank24
         rad(i) = 0.0d0
         eps(i) = 0.0d0
         rad4(i) = 0.0d0
         eps4(i) = 0.0d0
         reduct(i) = 0.0d0
         chg(i) = 0.0d0
         polr(i) = 0.0d0
         athl(i) = 0.0d0
         do j = 1, maxval
            pgrp(j,i) = 0
         end do
      end do
      do i = 1, maxclass
         do j = 1, 2
            stbn(j,i) = 0.0d0
         end do
         do j = 1, 3
            anan(j,i) = 0.0d0
         end do
         electron(i) = 0.0d0
         ionize(i) = 0.0d0
         repulse(i) = 0.0d0
      end do
      do i = 1, maxbio
         biotyp(i) = 0
      end do
c
c     initialize values of some MMFF-specific parameters
c
      do i = 1, mmffmaxnb0
         mmff_kb0(i) = 1000.0d0
         mmff_b0(i) = 1000.0d0
         mmffbs0(i) = blank8
      enddo
      do i = 1, mmffmaxnb1
         mmff_kb1(i) = 1000.0d0
         mmff_b1(i) = 1000.0d0
         mmffbs1(i) = blank8
      enddo
      do i = 1, mmffmaxnsb
         stbn_abc0(i) = 1000.0d0
         stbn_cba0(i) = 1000.0d0
         stbn_abc1(i) = 1000.0d0
         stbn_cba1(i) = 1000.0d0
         stbn_abc2(i) = 1000.0d0
         stbn_cba2(i) = 1000.0d0
         stbn_abc3(i) = 1000.0d0
         stbn_cba3(i) = 1000.0d0
         stbn_abc4(i) = 1000.0d0
         stbn_cba4(i) = 1000.0d0
         stbn_abc5(i) = 1000.0d0
         stbn_cba5(i) = 1000.0d0
         stbn_abc6(i) = 1000.0d0
         stbn_cba6(i) = 1000.0d0
         stbn_abc7(i) = 1000.0d0
         stbn_cba7(i) = 1000.0d0
         stbn_abc8(i) = 1000.0d0
         stbn_cba8(i) = 1000.0d0
         stbn_abc9(i) = 1000.0d0
         stbn_cba9(i) = 1000.0d0
         stbn_abc10(i) = 1000.0d0
         stbn_cba10(i) = 1000.0d0
         stbn_abc11(i) = 1000.0d0
         stbn_cba11(i) = 1000.0d0
      enddo 
      do i = 1, 100
         do j = 1, 100
            bci(j,i) = 1000.0d0
            bci_1(j,i) = 1000.0d0
         end do
      end do
      do i = 1, mmffmaxna
         mmff_ka0(i) = 1000.0d0
         mmff_ka1(i) = 1000.0d0
         mmff_ka2(i) = 1000.0d0
         mmff_ka3(i) = 1000.0d0
         mmff_ka4(i) = 1000.0d0
         mmff_ka5(i) = 1000.0d0
         mmff_ka6(i) = 1000.0d0
         mmff_ka7(i) = 1000.0d0
         mmff_ka8(i) = 1000.0d0
         mmff_ang0(i) = 1000.0d0
         mmff_ang1(i) = 1000.0d0
         mmff_ang2(i) = 1000.0d0
         mmff_ang3(i) = 1000.0d0
         mmff_ang4(i) = 1000.0d0
         mmff_ang5(i) = 1000.0d0
         mmff_ang6(i) = 1000.0d0
         mmff_ang7(i) = 1000.0d0
         mmff_ang8(i) = 1000.0d0
      end do
      do i = 1, maxele
         mmffanger(i,1) = 1000.0d0
         mmffanger(i,2) = 1000.0d0
      enddo
      do i = 1, maxnt
         kt(i) = blank16
         kt_1(i) = blank16
         kt_2(i) = blank16
         t1(1,i) = 1000.0d0
         t1(2,i) = 1000.0d0
         t2(1,i) = 1000.0d0
         t2(2,i) = 1000.0d0
         t3(1,i) = 1000.0d0
         t3(2,i) = 1000.0d0
         t1_1(1,i) = 1000.0d0
         t1_1(2,i) = 1000.0d0
         t2_1(1,i) = 1000.0d0
         t2_1(2,i) = 1000.0d0
         t3_1(1,i) = 1000.0d0
         t3_1(2,i) = 1000.0d0
         t1_2(1,i) = 1000.0d0
         t1_2(2,i) = 1000.0d0
         t2_2(1,i) = 1000.0d0
         t2_2(2,i) = 1000.0d0
         t3_2(1,i) = 1000.0d0
         t3_2(2,i) = 1000.0d0
      end do
      do i = 1, maxele
         mmfftorer(i,1) = 1000.0d0
         mmfftorer(i,2) = 1000.0d0
      end do
      do i = 1, 5
         do j = 1, maxtyp
            eqclass(j,i) = 1000
         end do
      end do
      do i = 1, 6
         do j = 1, maxtyp
            mmffarom(j,i) = 0
            mmffaromc(j,i) = 0
            mmffaroma(j,i) = 0
         end do
      end do
c
c     initialize values of Modified MMFF parameters
c
      do i = 1, maxnstrlf
         strlfid(i) = blank12
         do j = 1, 3
            strlfcon(i,j) = 0.0d0
         enddo
         strlfeq(i) = 0.0d0
      enddo
      do i = 1, maxnanglf
         anglfid(i) = blank16
         do j = 1, 3
            anglfcon(i,j) = 0.0d0
         enddo
         anglfeq(i) = 0.0d0
      enddo
      do i = 1, maxntorlf
         torlfid(i) = blank20
         do j = 1, 6
            torlfcon(i,j) = 0.0d0
            torlfph(i,j) = 0.0d0
         enddo
      enddo
c
c     initialize values of LFMM-related parameters
c
      mlrestlow = 0.3d0
      mlrestup = 1.2d0
      mlrestw8t = 1.0d8
      do i = 1, maxnhalf
         halfid(i) = blank8
         halfcon(i) = 0.0d0
         halfeq(i) = 0.0d0
      enddo
      do i = 1, maxnmolf
         molfid(i) = blank8
         molfcon(i,1) = 0.0d0
         molfcon(i,2) = 0.0d0
         molfeq(i) = 0.0d0
      enddo
      do i = 1, maxnlllf
         lllfid(i) = blank8
         lllfcon(i,1) = 0.0d0
         lllfcon(i,2) = 0.0d0
      enddo
      do i = 1, maxnvwlf
         vwlfid(i) = blank8
         vwlfcon(i,1) = 0.0d0
         vwlfcon(i,2) = 0.0d0
         vwlfcon(i,3) = 0.0d0
      enddo
      do i = 1, maxnaom
         aomsigid(i) = blank8
         aompixid(i) = blank8
         aompiyid(i) = blank8
         aomxdsid(i) = blank8
         epairid(i) = blank8
         do j = 1, 7
            aomsig(i,j) = 0.0d0
            aomxds(i,j) = 0.0d0
            aompix(i,j) = 0.0d0
            aompiy(i,j) = 0.0d0
            epairpar(i,j) = 0.0d0
         enddo
      enddo
c
c     set default control options for LFMM calculation
c
      do i = 1, maxlfmm
         metlfmm(i) = -1
         elelfmm(i) = -1
         spinlfmm(i) = -1
         sdlfsedone(i) = .false.
      enddo
      shajorlfse = 1
      barysjlfse = 2
      anlderlfse = 1
      stpzlfse = 0.02d0
      displfse = 2
      doderlfse = 0
      scnderlfse = 0
      stbenglfse = 0
      debuglfse = 0
      dosdfdout = .false.
      updtlfshess = .true.
      reproduceold = .false.
c
c     initialize values of classes to be treated by 12-10 vdw potential
c
      do i = 1, maxpsnb
         kisnb(1,i) = 0
         kisnb(2,i) = 0
      enddo
c
c     set default control parameters for local geometry terms
c
      bndtyp = 'HARMONIC'
      bndunit = 1.0d0
      cbnd = 0.0d0
      qbnd = 0.0d0
      angunit = 1.0d0 / radian**2
      cang = 0.0d0
      qang = 0.0d0
      pang = 0.0d0
      sang = 0.0d0
      stbnunit = 1.0d0 / radian
      ureyunit = 1.0d0
      cury = 0.0d0
      qury = 0.0d0
      aaunit = 1.0d0 / radian**2
      opbtyp = 'W-D-C'
      opbunit = 1.0d0 / radian**2
      copb = 0.0d0
      qopb = 0.0d0
      popb = 0.0d0
      sopb = 0.0d0
      opdunit = 1.0d0
      copd = 0.0d0
      qopd = 0.0d0
      popd = 0.0d0
      sopd = 0.0d0
      idihunit = 1.0d0
      itorunit = 1.0d0
      torsunit = 1.0d0
      ptorunit = 1.0d0
      storunit = 1.0d0
      ttorunit = 1.0d0
c
c     set default control parameters for van der Waals terms
c
      vdwindex = 'CLASS'
      vdwtyp = 'LENNARD-JONES'
      radrule = 'ARITHMETIC'
      radtyp = 'R-MIN'
      radsiz = 'RADIUS'
      epsrule = 'GEOMETRIC'
      gausstyp = 'NONE'
      ngauss = 0
      abuck = 0.0d0
      bbuck = 0.0d0
      cbuck = 0.0d0
      ghal = 0.12d0
      dhal = 0.07d0
      v2scale = 0.0d0
      v3scale = 0.0d0
      v4scale = 1.0d0
      v5scale = 1.0d0
      use_vcorr = .false.
c
c     set default control parameters for charge-charge terms
c
      electric = coulomb
      dielec = 1.0d0
      ebuffer = 0.0d0
      c2scale = 0.0d0
      c3scale = 0.0d0
      c4scale = 1.0d0
      c5scale = 1.0d0
      neutnbr = .false.
      neutcut = .false.
      rdepdielec = .false.
      tapertype = 1
c
c     set default control parameters for polarizable multipoles
c
      m2scale = 0.0d0
      m3scale = 0.0d0
      m4scale = 1.0d0
      m5scale = 1.0d0
      p2scale = 0.0d0
      p3scale = 0.0d0
      p4scale = 1.0d0
      p5scale = 1.0d0
      p41scale = 0.5d0
c
c     set default control parameters for induced dipoles
c
      poltyp = 'MUTUAL'
      politer = 500
      poleps = 0.000001d0
      udiag = 2.0d0
      d1scale = 0.0d0
      d2scale = 1.0d0
      d3scale = 1.0d0
      d4scale = 1.0d0
      u1scale = 1.0d0
      u2scale = 1.0d0
      u3scale = 1.0d0
      u4scale = 1.0d0
c
c     set default control parameters for implicit solvation
c
      solvtyp = blank8
      borntyp = blank8
c
c     set default control parameters for reaction field
c
      rfsize = 1000000.0d0
      rfbulkd = 80.0d0
      rfterms = 1
c
c     set default type of functional form for torsional restraints
c
      use_tfix2 = .false.
      return
      end
