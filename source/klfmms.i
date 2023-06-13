c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  klfmms.i  --  force field parameter storage for LFMM    ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxnstrlf maximum number of polynomial bond stretching parameters
c     maxnanglf maximum number of polynomial angle bending parameters
c     maxntorlf maximum number of polynomial torsion parameters
c     maxnhalf  maximum number of harmonic metal-ligand parameters
c     maxnmolf  maximum number of Morse metal-ligand parameters
c     maxnlllf  maximum number of ligand-ligand repulsion parameters
c     maxnvwlf  maximum number of van der Waals-like (Lennard-Jones-like) parameters
c     maxnaom   maximum number of angular overlap parameters
c
c     strlfid  string of atom classes and bond type for polynomial bond stretching term
c     strlfcon force constant parameters for polynomial bond stretching term
c     strlfeq  bond length parameters for polynomial bond stretching term
c 
c     anglfid  string of atom classes and angle type for polynomial angle bending term
c     anglfcon force constant parameters for polynomial angle bending term
c     anglfeq  bond angle for polynomial angle bending term
c
c     torlfid  string of atom classes and torsion type for polynomial torsion twist term
c     torlfcon force constant parameters for polynomial torsion twist term
c     torlfph  phase parameters for polynomial torsion twist term
c
c     halfid   string of atom classes for harmonic metal-ligand terms
c     halfcon  force constant parameters for harmonic metal-ligand terms
c     halfeq   bond length parameters for harmonic metal-ligand terms
c
c     molfid   string of atom classes for Morse metal-ligand terms
c     molfcon  parameters for Morse metal-ligand terms
c     molfeq   bond length parameters for Morse metal-ligand terms
c
c     mlrestlow lower distance limit of M-L restraint for Morse
c     mlrestup  upper distance limit of M-L restraint for Morse
c     mlrestw8t weight of M-L restraint for Morse
c
c     lllfid   string of atom classes for ligand-ligand repulsion terms
c     lllfcon  parameters for ligand-ligand repulsion terms
c
c     vwlfid   string of atom classes for van der Waals-like ligand repulsion
c     vwlfcon  force constant parameters for van der Waals-like ligand repulsion
c
c     aomsigid string of atom classes for angular overlap model - sigma terms
c     aomsig   parameters sigma for angular overlap terms
c     aompixid string of atom classes for angular overlap model - pix terms
c     aompix   parameters pi-x for angular overlap terms
c     aompiyid string of atom classes for angular overlap model - piy terms
c     aompiy   parameters pi-y for angular overlap terms
c     aomxdsid string of atom classes for angular overlap model - d-s mixing
c     aomxds   parameters for d-x mixing
c
c     epairid  string of atom classes for electron pairing energy terms
c     epairpar parameters for electron pairing energy terms
c
c
      integer maxnstrlf,maxnanglf,maxntorlf
      integer maxnhalf,maxnmolf,maxnlllf
      integer maxnvwlf,maxnaom
      parameter (maxnstrlf=4000)
      parameter (maxnanglf=4000)
      parameter (maxntorlf=4000)
      parameter (maxnhalf=100)
      parameter (maxnmolf=100)
      parameter (maxnlllf=100)
      parameter (maxnvwlf=100)
      parameter (maxnaom=100)
      real*8 strlfcon,strlfeq
      real*8 anglfcon,anglfeq
      real*8 torlfcon,torlfph
      real*8 halfcon,halfeq
      real*8 molfcon,molfeq
      real*8 lllfcon
      real*8 vwlfcon
      real*8 aomsig,aompix,aompiy,aomxds
      real*8 epairpar
      real*8 mlrestlow,mlrestup,mlrestw8t
      character*8 halfid,molfid,lllfid
      character*8 vwlfid,epairid
      character*8 aomsigid,aompixid
      character*8 aompiyid,aomxdsid
      character*12 strlfid
      character*16 anglfid
      character*20 torlfid
      common /klfmms/ mlrestlow,mlrestup,mlrestw8t,
     &                strlfid(maxnstrlf),strlfcon(maxnstrlf,3),
     &                strlfeq(maxnstrlf),
     &                anglfid(maxnanglf),anglfcon(maxnanglf,3),
     &                anglfeq(maxnanglf),
     &                torlfid(maxntorlf),torlfcon(maxntorlf,6),
     &                torlfph(maxntorlf,6),
     &                halfid(maxnhalf),halfcon(maxnhalf),
     &                halfeq(maxnhalf),
     &                molfid(maxnmolf),molfcon(maxnmolf,2),
     &                molfeq(maxnmolf),
     &                lllfid(maxnlllf),lllfcon(maxnlllf,2),
     &                vwlfid(maxnvwlf),vwlfcon(maxnvwlf,3),
     &                aomsigid(maxnaom),aompixid(maxnaom),
     &                aompiyid(maxnaom),aomxdsid(maxnaom),
     &                aomsig(maxnaom,7),aomxds(maxnaom,7),
     &                aompix(maxnaom,7),aompiy(maxnaom,7),
     &                epairid(maxnaom),epairpar(maxnaom,7)
