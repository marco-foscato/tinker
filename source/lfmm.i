c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  lfmm.i  --  parameters for LFMM calculation in the  ##
c     ##              current structure                       ##
c     ##                                                      ##
c     ##########################################################
c
c 
c     nhalf    total number of harmonic terms from LFMM
c     ihalf    numbers of the atoms in each harmonic terms from LFMM
c     lhalf    ideal bond length values in Angstroms for each harmonic terms from LFMM
c     khalf    force constants for each harmonic terms from LFMM
c
c     nmolf    total number of Morse terms from LFMM
c     imolf    numbers of the atoms in each Morse terms from LFMM
c     lmolf    ideal bond length values in Angstroms for each Morse terms from LFMM
c     bmolf    parameter 'b' for each Morse potential from LFMM
c     amolf    parameter 'a' for each Morse potential from LFMM
c
c     lowmolf  lower distance limit of M-L restraint for Morse
c     upmolf   upper distance limit of M-L restraint for Morse
c     w8tmolf  weight of M-L restraint for Morse
c
c     nlllf    total number of ligand-ligand interactions from LFMM
c     illlf    numbers of the atoms in each ligand-ligand interactions from LFMM
c     alllf    parameter 'a' for each ligand-ligand interactions from LFMM
c     mlllf    parameter 'm' for each ligand-ligand interactions from LFMM
c
c     nvwlf    total number of Lennard-Jones-like terms from LFMM
c     ivwlf    numbers of the atoms in each Lennard-Jones-like term from LFMM
c     mvwlf    exponent for each Lennard-Jones-like terms from LFMM
c     avwlf    parameter 'b' for each Lennard-Jones-like terms from LFMM
c     bvwlf    parameter 'a' for each Lennard-Jones-like terms from LFMM
c
c     npair    total number of electron pairing terms
c     ipair    numbers of the atoms in each electron pairing term
c     ppair    parameters for electron pair repulsion
c     wpair    pre-factor for electron pair repulsion
c
c     nlfse    total number of ligand field stabilization energy terms
c     imlfse   numbers of the central atom in each ligand field term
c     illfse   numbers of the coordinating atoms in each ligand field term
c     iss1     numbers of the atoms in the first group of second shell atoms for each ligand in each LFMM center 
c     iss2     numbers of the atoms in the second group of second shell atoms for each ligand in each LFMM center
c     nligs    total number of ligands in each LFMM center
c     nss1     number of second shell, non-metal atoms in the 1st group neighbours for each ligand in each LFMM center
c     nss2     number of second shell, non-metal atoms in the 2nd group neighbours for each ligand in each LFMM center
c     totaom   total number of atoms in each ligand field term
c     aom2tnk  reordered list of atoms; aom2tnk(idcenter,i) = j -> j-th atom in TINKER list is the i-th in AOM region
c     ligunqtyps list of unique ligand types with a corresponding ligand for each LFMM center
c     ligtyp   indexes of ligand type in ligtypes for each ligand in each LFMM center
c     nelelfse number of d-electrons per each ligand field stabilization energy term
c     smlfse   spin state flag (0=low; 1=high; 2=intermediate) per each lfse
c
c     esig     AOM parameters sigma contribution
c     epix     AOM parameters pi-x contribution
c     epiy     AOM parameters pi-y contribution
c     exds     AOM parameters d/s mixing contribution
c
c
      integer nhalf,nmolf,nlllf
      integer nvwlf,mvwlf
      integer npair,nlfse,nligs
      integer ihalf,imolf,illlf
      integer ivwlf,ipair
      integer imlfse,illfse
      integer iss1,iss2
      integer nss1,nss2
      integer aom2tnk,totaom
      integer ligtyp,ligunqtyps
      integer nelelfse,smlfse
      real*8 lhalf,khalf
      real*8 lmolf,bmolf,amolf
      real*8 lowmolf,upmolf,w8tmolf
      real*8 alllf,mlllf
      real*8 avwlf,bvwlf
      real*8 esig,epix,epiy,exds
      real*8 ppair,wpair

      common /lfmm/ lhalf(maxlfmmtrm),khalf(maxlfmmtrm),
     &              lmolf(maxlfmmtrm),bmolf(maxlfmmtrm),
     &              lowmolf,upmolf,w8tmolf,
     &              amolf(maxlfmmtrm),
     &              alllf(maxlfmmtrm),mlllf(maxlfmmtrm),
     &              avwlf(maxlfmmtrm),bvwlf(maxlfmmtrm),
     &              esig(maxlfmm,maxlig,7),epix(maxlfmm,maxlig,7),
     &              epiy(maxlfmm,maxlig,7),exds(maxlfmm,maxlig,7),
     &              ppair(maxlfmmtrm,7),wpair(maxlfmmtrm),
     &              aom2tnk(maxlfmm,maxcrd),
     &              totaom(maxlfmm),
     &              nmolf,nhalf,nlllf,nvwlf,npair,nlfse,
     &              nligs(maxlfmm),
     &              nss1(maxlfmm,maxlig),nss2(maxlfmm,maxlig),
     &              iss1(maxlfmm,maxlig),
     &              iss2(maxlfmm,maxlig),
     &              ihalf(2,maxlfmmtrm),imolf(2,maxlfmmtrm),
     &              illlf(2,maxlfmmtrm),ivwlf(2,maxlfmmtrm),
     &              ipair(2,maxlfmmtrm),
     &              imlfse(maxlfmm),illfse(maxlfmm,maxlig),
     &              mvwlf(maxlfmmtrm),
     &              nelelfse(maxlfmm),smlfse(maxlfmm),
     &              ligtyp(maxlfmm,maxlig),ligunqtyps(maxlfmm,maxlig,2)
