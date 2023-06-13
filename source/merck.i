c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2007 by Nicolas Staelens & Jay W. Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  merck.i  --  parameters specific for MMFF force field  ##
c     ##                                                         ##
c     #############################################################
c
c
c     mmffmaxnb0 maximum number of bond stretch parameter entries 
c                having MMFF Bond Type 0
c     mmffmaxnb1 maximum number of bond stretch parameter entries
c                having MMFF Bond Type 1
c     mmffmaxbt1 maximum number of bonds having MMFF Bond Type 1 
c
c     bt_1       atom pairs having MMFF Bond Type 1
c     nlignes    number of atom pairs having MMFF Bond Type 1
c     eqclass    table of atom class equivalencies used to find
c                default parameters if explicit values are missing
c                (see J. Comput. Chem., 17, 490-519, 95, Table IV)
c     crd        number of attached neighbors    |
c     val        valency value                   |  see T. A. Halgren,
c     pilp       if 0, no lone pair              |  J. Comput. Chem.,
c                if 1, one or more lone pair(s)  |  17, 616-645 (1995)
c     mltb       multibond indicator             |
c     arom       aromaticity indicator           |
c     lin        linearity indicator             |
c     sbmb       single- vs multiple-bond flag   |
c     mmffarom   aromatic rings parameters
c     mmffaromc  cationic aromatic rings parameters
c     mmffaroma  anionic aromatic rings parameters
c     mmffbs0    string of atom classes for bond stretch Type 0
c     mmffbs1    string of atom classes for bond stretch Type 1
c
c
      integer mmffmaxnb0
      integer mmffmaxnb1
      integer mmffmaxbt1
      integer bt_1
      integer nlignes
      integer eqclass
      integer crd,val,pilp,mltb
      integer arom,lin,sbmb
      integer mmffarom
      integer mmffaromc
      integer mmffaroma
      parameter (mmffmaxnb0=2000)
      parameter (mmffmaxnb1=2000)
      parameter (mmffmaxbt1=500)
      character*8 mmffbs0,mmffbs1
      common /merck1/ bt_1(mmffmaxbt1,2),nlignes,
     &                eqclass(maxtyp,5),crd(maxclass),
     &                val(maxclass),pilp(maxclass),
     &                mltb(maxclass),arom(maxclass),
     &                lin(maxclass),sbmb(maxclass),
     &                mmffarom(maxtyp,6),
     &                mmffaromc(maxtyp,6),mmffaroma(maxtyp,6),
     &                mmffbs0(mmffmaxnb0),mmffbs1(mmffmaxnb1)
c
c
c     mmff_kb0  bond force constant for bond stretch with Bond Type 0
c     mmff_kb1  bond force constant for bond stretch with Bond Type 1
c     mmff_b0   bond length value for bond stretch with Bond Type 0
c     mmff_b1   bond length value for bond stretch with Bond Type 1
c     rad0      covalent atomic radius for empirical bond rules
c     paulel    Pauling electronegativities for empirical bond rules
c     r0ref     reference bond length for empirical bond rules
c     kbref     reference force constant for empirical bond rules
c
c
      real*8 mmff_kb0,mmff_kb1
      real*8 mmff_b0,mmff_b1
      real*8 rad0,r0ref
      real*8 kbref,paulel
      common /merck2/ mmff_kb0(mmffmaxnb0),
     &                mmff_kb1(mmffmaxnb1),
     &                mmff_b0(mmffmaxnb0),
     &                mmff_b1(mmffmaxnb1),
     &                rad0(maxele),paulel(maxele),
     &                r0ref(maxele,maxele),
     &                kbref(maxele,maxele)
c
c
c     mmffmaxna  maximum number of angle bending parameters
c
c     mmffangcl0 string of atom classes for angle bend
c     mmffangcl1 string of atom classes for angle bend with one bond of Type 1
c     mmffangcl2 string of atom classes for angle bend with both bonds of Type 1
c     mmffangcl3 string of atom classes for angle bend for 3-membered ring
c     mmffangcl4 string of atom classes for angle bend for 4-membered ring
c     mmffangcl5 string of atom classes for angle bend for 3-ring and one Bond Type 1
c     mmffangcl6 string of atom classes for angle bend for 3-ring and both Bond Type 1
c     mmffangcl7 string of atom classes for angle bend for 4-ring and one Bond Type 1
c     mmffangcl8 string of atom classes for angle bend for 4-ring and both Bond Type 1
c
c     mmff_ka0    angle force constant for triples of atom classes
c     mmff_ka1    angle force constant with one bond of Type 1
c     mmff_ka2    angle force constant with both bonds of Type 1
c     mmff_ka3    angle force constant for 3-membered ring
c     mmff_ka4    angle force constant for 4-membered ring
c     mmff_ka5    angle force constant for 3-ring and one Bond Type 1
c     mmff_ka6    angle force constant for 3-ring and both Bond Type 1
c     mmff_ka7    angle force constant for 4-ring and one Bond Type 1
c     mmff_ka8    angle force constant for 4-ring and both Bond Type 1
c     mmff_ang0   ideal bond angle for triples of atom classes
c     mmff_ang1   ideal bond angle with one bond of Type 1
c     mmff_ang2   ideal bond angle with both bonds of Type 1
c     mmff_ang3   ideal bond angle for 3-membered ring
c     mmff_ang4   ideal bond angle for 4-membered ring
c     mmff_ang5   ideal bond angle for 3-ring and one Bond Type 1
c     mmff_ang6   ideal bond angle for 3-ring and both Bond Type 1
c     mmff_ang7   ideal bond angle for 4-ring and one Bond Type 1
c     mmff_ang8   ideal bond angle for 4-ring and both Bond Type 1
c
c     mmffanger   parameters for empirical angle rule 
c
      integer mmffmaxna
      parameter (mmffmaxna=2000)
      real*8 mmff_ka0,mmff_ka1,mmff_ka2
      real*8 mmff_ka3,mmff_ka4,mmff_ka5
      real*8 mmff_ka6,mmff_ka7,mmff_ka8
      real*8 mmff_ang0,mmff_ang1,mmff_ang2
      real*8 mmff_ang3,mmff_ang4,mmff_ang5
      real*8 mmff_ang6,mmff_ang7,mmff_ang8
      real*8 mmffanger
      character*12 mmffangcl0,mmffangcl1,mmffangcl2
      character*12 mmffangcl3,mmffangcl4,mmffangcl5
      character*12 mmffangcl6,mmffangcl7,mmffangcl8
      common /merck3/  mmff_ka0(mmffmaxna),mmff_ka1(mmffmaxna),
     &                 mmff_ka2(mmffmaxna),mmff_ka3(mmffmaxna),
     &                 mmff_ka4(mmffmaxna),mmff_ka5(mmffmaxna),
     &                 mmff_ka6(mmffmaxna),mmff_ka7(mmffmaxna),
     &                 mmff_ka8(mmffmaxna),
     &                 mmff_ang0(mmffmaxna),mmff_ang1(mmffmaxna),
     &                 mmff_ang2(mmffmaxna),mmff_ang3(mmffmaxna),
     &                 mmff_ang4(mmffmaxna),mmff_ang5(mmffmaxna),
     &                 mmff_ang6(mmffmaxna),mmff_ang7(mmffmaxna),
     &                 mmff_ang8(mmffmaxna),
     &                 mmffangcl0(mmffmaxna),mmffangcl1(mmffmaxna),
     &                 mmffangcl2(mmffmaxna),mmffangcl3(mmffmaxna),
     &                 mmffangcl4(mmffmaxna),mmffangcl5(mmffmaxna),
     &                 mmffangcl6(mmffmaxna),mmffangcl7(mmffmaxna),
     &                 mmffangcl8(mmffmaxna),
     &                 mmffanger(maxele,2)
c
c
c     mmffmaxnsb  maximum number of stretch bend parameters
c
c     Stretch-Bend Type 0
c     mmffsbabclc0 string of atom classes for stretch bend parameters
c     mmffscbaclc0 string of atom classes for stretch bend parameters
c     stbn_abc0    stretch-bend parameters for A-B-C atom classes
c     stbn_cba0    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 1  (A-B is Bond Type 1)
c     mmffsbabclc1 string of atom classes for stretch bend parameters
c     mmffscbaclc1 string of atom classes for stretch bend parameters
c     stbn_abc1    stretch-bend parameters for A-B-C atom classes
c     stbn_cba1    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 2  (B-C is Bond Type 1) 
c     mmffsbabclc2 string of atom classes for stretch bend parameters
c     mmffscbaclc2 string of atom classes for stretch bend parameters
c     stbn_abc2    stretch-bend parameters for A-B-C atom classes
c     stbn_cba2    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 3  (A-B and B-C are Bond Type 1) 
c     mmffsbabclc3 string of atom classes for stretch bend parameters
c     mmffscbaclc3 string of atom classes for stretch bend parameters
c     stbn_abc3    stretch-bend parameters for A-B-C atom classes
c     stbn_cba3    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 4  (both Bond Types 0, 4-membered ring)
c     mmffsbabclc4 string of atom classes for stretch bend parameters
c     mmffscbaclc4 string of atom classes for stretch bend parameters
c     stbn_abc4    stretch-bend parameters for A-B-C atom classes
c     stbn_cba4    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 5  (both Bond Types 0, 3-membered ring)
c     mmffsbabclc5 string of atom classes for stretch bend parameters
c     mmffscbaclc5 string of atom classes for stretch bend parameters
c     stbn_abc5    stretch-bend parameters for A-B-C atom classes
c     stbn_cba5    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 6  (A-B is Bond Type 1, 3-membered ring)
c     mmffsbabclc6 string of atom classes for stretch bend parameters
c     mmffscbaclc6 string of atom classes for stretch bend parameters
c     stbn_abc6    stretch-bend parameters for A-B-C atom classes
c     stbn_cba6    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 7  (B-C is Bond Type 1, 3-membered ring)
c     mmffsbabclc7 string of atom classes for stretch bend parameters
c     mmffscbaclc7 string of atom classes for stretch bend parameters
c     stbn_abc7    stretch-bend parameters for A-B-C atom classes
c     stbn_cba7    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 8  (both Bond Types 1, 3-membered ring)
c     mmffsbabclc8 string of atom classes for stretch bend parameters
c     mmffscbaclc8 string of atom classes for stretch bend parameters
c     stbn_abc8    stretch-bend parameters for A-B-C atom classes
c     stbn_cba8    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 9  (A-B is Bond Type 1, 4-membered ring)
c     mmffsbabclc9 string of atom classes for stretch bend parameters
c     mmffscbaclc9 string of atom classes for stretch bend parameters
c     stbn_abc9    stretch-bend parameters for A-B-C atom classes
c     stbn_cba9    stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 10  (B-C is Bond Type 1, 4-membered ring)
c     mmffsbabclc10 string of atom classes for stretch bend parameters
c     mmffscbaclc10 string of atom classes for stretch bend parameters
c     stbn_abc10   stretch-bend parameters for A-B-C atom classes
c     stbn_cba10   stretch-bend parameters for C-B-A atom classes
c
c     Stretch-Bend Type 11  (both Bond Types 1, 4-membered ring)
c     mmffsbabclc11 string of atom classes for stretch bend parameters
c     mmffscbaclc11 string of atom classes for stretch bend parameters
c     stbn_abc11   stretch-bend parameters for A-B-C atom classes
c     stbn_cba11   stretch-bend parameters for C-B-A atom classes
c
c     defstbn_abc  default stretch-bend parameters for A-B-C classes
c     defstbn_cba  default stretch-bend parameters for C-B-A classes
c
c
      integer mmffmaxnsb
      parameter (mmffmaxnsb=2000)
      real*8 stbn_abc0,stbn_cba0
      real*8 stbn_abc1,stbn_cba1
      real*8 stbn_abc2,stbn_cba2
      real*8 stbn_abc3,stbn_cba3
      real*8 stbn_abc4,stbn_cba4
      real*8 stbn_abc5,stbn_cba5
      real*8 stbn_abc6,stbn_cba6
      real*8 stbn_abc7,stbn_cba7
      real*8 stbn_abc8,stbn_cba8
      real*8 stbn_abc9,stbn_cba9
      real*8 stbn_abc10,stbn_cba10
      real*8 stbn_abc11,stbn_cba11
      real*8 defstbn_abc,defstbn_cba
      character*12 mmffsbabclc0,mmffscbaclc0
      character*12 mmffsbabclc1,mmffscbaclc1
      character*12 mmffsbabclc2,mmffscbaclc2
      character*12 mmffsbabclc3,mmffscbaclc3
      character*12 mmffsbabclc4,mmffscbaclc4
      character*12 mmffsbabclc5,mmffscbaclc5
      character*12 mmffsbabclc6,mmffscbaclc6
      character*12 mmffsbabclc7,mmffscbaclc7
      character*12 mmffsbabclc8,mmffscbaclc8
      character*12 mmffsbabclc9,mmffscbaclc9
      character*12 mmffsbabclc10,mmffscbaclc10
      character*12 mmffsbabclc11,mmffscbaclc11

      common /merck4/ stbn_abc0(mmffmaxnsb),stbn_cba0(mmffmaxnsb),
     &                stbn_abc1(mmffmaxnsb),stbn_cba1(mmffmaxnsb),
     &                stbn_abc2(mmffmaxnsb),stbn_cba2(mmffmaxnsb),
     &                stbn_abc3(mmffmaxnsb),stbn_cba3(mmffmaxnsb),
     &                stbn_abc4(mmffmaxnsb),stbn_cba4(mmffmaxnsb),
     &                stbn_abc5(mmffmaxnsb),stbn_cba5(mmffmaxnsb),
     &                stbn_abc6(mmffmaxnsb),stbn_cba6(mmffmaxnsb),
     &                stbn_abc7(mmffmaxnsb),stbn_cba7(mmffmaxnsb),
     &                stbn_abc8(mmffmaxnsb),stbn_cba8(mmffmaxnsb),
     &                stbn_abc9(mmffmaxnsb),stbn_cba9(mmffmaxnsb),
     &                stbn_abc10(mmffmaxnsb),stbn_cba10(mmffmaxnsb),
     &                stbn_abc11(mmffmaxnsb),stbn_cba11(mmffmaxnsb),
     &                defstbn_abc(0:4,0:4,0:4),
     &                defstbn_cba(0:4,0:4,0:4),
     &                mmffsbabclc0(mmffmaxnsb),
     &                mmffscbaclc0(mmffmaxnsb),
     &                mmffsbabclc1(mmffmaxnsb),
     &                mmffscbaclc1(mmffmaxnsb),
     &                mmffsbabclc2(mmffmaxnsb),
     &                mmffscbaclc2(mmffmaxnsb),
     &                mmffsbabclc3(mmffmaxnsb),
     &                mmffscbaclc3(mmffmaxnsb),
     &                mmffsbabclc4(mmffmaxnsb),
     &                mmffscbaclc4(mmffmaxnsb),
     &                mmffsbabclc5(mmffmaxnsb),
     &                mmffscbaclc5(mmffmaxnsb),
     &                mmffsbabclc6(mmffmaxnsb),
     &                mmffscbaclc6(mmffmaxnsb),
     &                mmffsbabclc7(mmffmaxnsb),
     &                mmffscbaclc7(mmffmaxnsb),
     &                mmffsbabclc8(mmffmaxnsb),
     &                mmffscbaclc8(mmffmaxnsb),
     &                mmffsbabclc9(mmffmaxnsb),
     &                mmffscbaclc9(mmffmaxnsb),
     &                mmffsbabclc10(mmffmaxnsb),
     &                mmffscbaclc10(mmffmaxnsb),
     &                mmffsbabclc11(mmffmaxnsb),
     &                mmffscbaclc11(mmffmaxnsb)
c
c
c     t1_1     torsional parameters for 1-fold, MMFF Torsion Type 1
c     t1_2     torsional parameters for 1-fold, MMFF Torsion Type 2
c     t2_1     torsional parameters for 2-fold, MMFF Torsion Type 1
c     t2_2     torsional parameters for 2-fold, MMFF Torsion Type 2
c     t3_1     torsional parameters for 3-fold, MMFF Torsion Type 1
c     t3_2     torsional parameters for 3-fold, MMFF Torsion Type 2
c     kt_1     string of classes for torsions, MMFF Torsion Type 1
c     kt_2     string of classes for torsions, MMFF Torsion Type 2
c
c     mmfftorer parameters for empirical torsion rules
c
      real*8 t1_1,t2_1,t3_1
      real*8 t1_2,t2_2,t3_2
      real*8 mmfftorer
      character*16 kt_1,kt_2
      common /merck5/ t1_1(2,0:maxtors),t2_1(2,0:maxtors),
     &                t3_1(2,0:maxtors),t1_2(2,0:maxtors),
     &                t2_2(2,0:maxtors),t3_2(2,0:maxtors),
     &                kt_1(0:maxtors),kt_2(0:maxtors),
     &                mmfftorer(maxele,2)
c
c
c     g        scale factors for calculation of MMFF eps
c     alph     atomic polarizabilities for calculation of MMFF eps
c     nn       effective number of valence electrons for MMFF eps
c     da       donor/acceptor atom classes
c
c
      real*8 g,alph,nn
      character*1 da
      common /merck6/ g(maxclass),alph(maxclass),nn(maxclass),
     &                da(maxclass)
c
c
c     bci      bond charge increments for building atom charges
c     bci_1    bond charge increments for MMFF Bond Type 1   
c     pbci     partial BCI for building missing BCIs
c     fcadj    formal charge adjustment factor
c
c
      real*8 bci,bci_1
      real*8 pbci,fcadj
      common /merck7/ bci(maxclass,maxclass),
     &                bci_1(maxclass,maxclass),
     &                pbci(maxclass),fcadj(maxclass)
