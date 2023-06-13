c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine klfmm  --  LFMM parameters assignment   ##
c     ##                                                     ##
c     #########################################################
c
c
c     "klfmm" assigns Ligand Field Molecular Mechanic (LFMM) 
c     parameters to transition metal center as specified in the key file.
c
c
      subroutine klfmm
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'bond.i'
      include 'angle.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'klfmms.i'
      include 'lfmmset.i'
      include 'lfmm.i'
      include 'potent.i'
      integer i,j,k,l,nl
      integer ia,ib,ic,im
      integer ita,itb,itc,itm
      integer next,size
      integer mvw_a,mvw_c
      real*8 wtmp
      real*8 all_a,all_c
      real*8 mll_a,mll_c
      real*8 avw_a,avw_c
      real*8 bve_a,bve_c
      character*4 pa,pb,pc,pm
      character*8 blank,pt
      character*20 keyword
      character*120 record,string,blank120
      logical header
      logical foundm,foundh
      logical foundll_a,foundvw_a
      logical foundll_c,foundvw_c
      logical donell,donevw
      logical foundsig,foundpix
      logical foundpiy,foundxds
      logical foundpair
c
c
c     initializations
c
      blank = '        '
      blank120 = repeat (' ',120)
      nlfmm = 0
      nmolf = 0
      nhalf = 0
      nlllf = 0
      nvwlf = 0
      nlfse = 0
      npair = 0
      do i = 1, maxlfmm
         do j = 1, maxlig
            ligtyp(i,j) = -1
            ligunqtyps(i,j,1) = -1
            ligunqtyps(i,j,2) = 0
         enddo
      enddo
      lowmolf = mlrestlow
      upmolf = mlrestup
      w8tmolf = mlrestw8t
c
c     process LFMM keywords
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'LFMMCENTER ') then
            nlfmm = nlfmm + 1
            if (nlfmm .ge. maxlfmm) then
               write (iout,10)
   10          format (/,' KLFMM  --  Too many LFMM centers;',
     &                   ' Increase MAXLFMM')
               call fatal
            endif
            read (string,*,err=20,end=20) ia,ib,ic
   20       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Defining LFMM centers: ',
     &                    //,5x,'Atom Number',6x,'d-Elect.',6x,'Spin',/)
               end if
               write(iout,40) ia,ib,ic
   40          format (8x,i4,15x,i2,11x,i1)
            end if
c
c           Check atom index
c
            if (ia .lt. 1) then
               write (iout,50) 'metal ID'
   50          format (/,' KLFMM  --  Incorrect definition of',
     &                   ' a LFMM center; Check input (',a,')')
               call fatal
c
c           Check number of electrons
c
            else if ((ib .lt. 0) .or. (ib .gt. 10)) then
               write (iout,50) 'no. electrons'
               call fatal
c
c           Check spin state indicator
c
            else if ((ic .lt. 0) .or. (ic .gt. 2)) then
               write (iout,50) 'spin state'
               call fatal
            end if
            metlfmm(nlfmm) = ia
            elelfmm(nlfmm) = ib
            spinlfmm(nlfmm) = ic
c
c        disable Shaffer and Jorgensen d/s-mixing
c
         else if (keyword(1:11) .eq. 'NO-SJ-LFSE ') then
            shajorlfse = 0
c
c        set baricenter of Shaffer and Jorgensen d/s-mixing
c
         else if (keyword(1:13) .eq. 'BARY-SJ-LFSE ') then
            read (string,*,err=69,end=69)  barysjlfse
c
c        set calculation of energy derivatives via finite differences
c
         else if (keyword(1:14) .eq. 'DERIV-FD-LFSE ') then
            anlderlfse = 0
c
c        set step size for finite difference
c
         else if (keyword(1:13) .eq. 'STEP-FD-LFSE ') then
            read (string,*,err=69,end=69)  stpzlfse
c
c        set displacement for finite difference
c
         else if (keyword(1:14) .eq. 'DISPL-FD-LFSE ') then
            read (string,*,err=69,end=69)  displfse
c        
c        ask Tinker to do finite difference for lfse 2nd derivatives
c
         else if (keyword(1:14) .eq. 'TNK-FDSD-LFSE ') then
            dosdfdout = .true.
c        
c        set use of stabilizzation energy formula
c
         else if (keyword(1:18) .eq. 'STB-ENG-FORM-LFSE ') then
            stbenglfse = 1
c
c        set flag for reproducing pre-2014 behaviour
c
         else if (keyword(1:19) .eq. 'REPRODUCE-PRE-2014 ') then
            reproduceold = .true.
c
c        set lower distance limit of M-L restraint for Morse
c
         else if (keyword(1:17) .eq. 'LOW-ML-RESTRAINT ') then
            read (string,*,err=69,end=69) lowmolf
c
c        set upper distance limit of M-L restraint for Morse
c
         else if (keyword(1:16) .eq. 'UP-ML-RESTRAINT ') then
            read (string,*,err=69,end=69) upmolf
c
c        set weight of M-L restraint for Morse
c
         else if (keyword(1:15) .eq. 'W-ML-RESTRAINT ') then
            read (string,*,err=69,end=69) w8tmolf
c
c        set debug printing of LFSE routines
c
         else if (keyword(1:11) .eq. 'DEBUG-LFSE ') then
            debuglfse = 1
         endif
  69     continue
      enddo
c
c     check identification of metal centers
c
      if (nlfmm .eq. 0) then
         write(iout,71)
  71     format (/,'  KLFMM  --  No metal center defined; '
     &             'Define values of LFMMCENTER in KEY file')
         call fatal
      end if
c
c     set bond contribution to zero if atom is LFMM center
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         do j = 1, nlfmm
            if ((ia .eq. metlfmm(j)) .or. (ib .eq. metlfmm(j))) then
               bk(i) = 0.0d0
               do k = 1, 3
                  bkpoly(i,k) = 0.0d0
               enddo               
            end if
         enddo
      enddo
c
c     set angle contributions to zero if central atom is LFMM center
c     
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         do j = 1, nlfmm
            if (ib .eq. metlfmm(j)) then
               ak(i) = 0.0d0
               anat(i) = 0.0d0
               afld(i) = 0.0d0
               do k = 1, 3
                  akpoly(i,k) = 0.0d0
               enddo
            end if
         enddo
      enddo 
c
c     assign parameters to metal-ligand bond terms
c
      if (use_lfmmhar .or. use_lfmmmor) then
         do i = 1, nlfmm
            ia = metlfmm(i)
            do j = 1, n12(ia)
               ib = i12(j,ia)
               ita = class(ia)
               itb = class(ib)
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               if (ita .le. itb) then
                  pt = pa//pb
               else
                  pt = pb//pa
               end if
c
c           Look for parameters for Morse potential
c
               foundm = .false.
               do k = 1, maxnmolf
                  if (molfid(k) .eq. pt) then
                     nmolf = nmolf + 1
                     imolf(1,nmolf) = ia
                     imolf(2,nmolf) = ib
                     lmolf(nmolf) = molfeq(k)
                     amolf(nmolf) = molfcon(k,2)
                     bmolf(nmolf) = molfcon(k,1)
                     foundm = .true.
                     exit
                  endif
               enddo
c
c           Look for parameters for Morse potential
c
               foundh = .false.
               do k = 1, maxnhalf
                  if (halfid(k) .eq. pt) then
                     nhalf = nhalf + 1
                     ihalf(1,nhalf) = ia
                     ihalf(2,nhalf) = ib
                     lhalf(nhalf) = halfeq(k)
                     khalf(nhalf) = halfcon(k)
                     foundh = .true.
                     exit
                  endif
               enddo
               if (.not. foundm .and. .not. foundh) then
                  write(iout,80) ita,itb
  80              format (/,' KLFMM  --  Unable to find LFMM parameters'
     &                    ' for metal-ligand',/,' bond stretching'
     &                    ' (classes ',i5,' and ',i5,');',/,
     &                    ' Check input and force field file.')
                  call fatal
               endif
            enddo
         enddo
      endif
c
c     assign parameters to ligand-ligand repulsion terms
c
      if (use_lfmmllr .or. use_lfmmvwl) then
         do i = 1, nlfmm
            im = metlfmm(i)
            itm = class(im)
            size = 4
            call numeral (itm,pm,size)
c
c        Search all as 1-3 ligand-ligand relations
c        
            do j = 1, n12(im)
               ia = i12(j,im)
               ita = class(ia)
               size = 4
               call numeral (ita,pa,size)
               if (ita .le. itm) then
                  pt = pa//pm
               else
                  pt = pm//pa
               end if
c
c           get preliminar parameters for first atom in 1-3 relation
c
               all_a = 0.0d0
               mll_a = 0.0d0
               foundll_a = .false.
               do l = 1, maxnlllf
                  if (lllfid(l) .eq. pt) then
                     all_a = lllfcon(l,1)
                     mll_a = lllfcon(l,2)
                     foundll_a = .true.
                     exit
                  endif
               enddo
               mvw_a = 0
               avw_a = 0.0d0
               bve_a = 0.0d0
               foundvw_a = .false.
               do l = 1, maxnvwlf
                  if (vwlfid(l) .eq. pt) then
                     mvw_a = int(vwlfcon(l,1))
                     avw_a = vwlfcon(l,2)
                     bve_a = vwlfcon(l,3)
                     foundvw_a = .true.
                     exit
                  endif
               enddo
               if (.not. foundll_a .and. .not. foundvw_a) then
                  write(iout,90) itm,ita
  90              format (/,' KLFMM  --  Unable to find LFMM parameters'
     &                    ' for ligand-ligand',/,' repultion term'
     &                    ' (classes ',i5,' and ',i5,');',/,
     &                    ' Check input and force field file.')
                  call fatal
               endif
c
c           identify second atom in 1-3 relation
c
               do k = j+1, n12(im)
                  ic = i12(k,im)
                  itc = class(ic)
                  size = 4
                  call numeral (itc,pc,size)
                  if (itc .le. itm) then
                     pt = pc//pm
                  else
                     pt = pm//pc
                  end if
c
c              get parameters second atom in 1-3 relation
c
                  all_c = 0.0d0
                  mll_c = 0.0d0
                  foundll_c = .false.
                  do l = 1, maxnlllf
                     if (lllfid(l) .eq. pt) then
                        all_c = lllfcon(l,1)
                        mll_c = lllfcon(l,2)
                        foundll_c = .true.
                        exit
                     endif
                  enddo
                  mvw_c = 0
                  avw_c = 0.0d0
                  bve_c = 0.0d0
                  foundvw_c = .false.
                  do l = 1, maxnvwlf
                     if (vwlfid(l) .eq. pt) then
                        mvw_c = int(vwlfcon(l,1))
                        avw_c = vwlfcon(l,2)
                        bve_c = vwlfcon(l,3)
                        foundvw_c = .true.
                        exit
                     endif
                  enddo
                  if (.not. foundll_c .and. .not. foundvw_c) then
                     write(iout,90) itm,itc
                     call fatal
                  endif
c
c              calculate combined parameters with mixing rule
c
                  if (foundll_a .and. foundll_c) then
                     nlllf = nlllf + 1
                     illlf(1,nlllf) = ia
                     illlf(2,nlllf) = ic
c
c                 geometric combination rule
c
                     alllf(nlllf) = SQRT(all_a * all_c)
                     mlllf(nlllf) = SQRT(mll_a * mll_c)
                     donell = .true.
                  endif
                  if (foundvw_a .and. foundvw_c) then
                     nvwlf = nvwlf + 1
                     ivwlf(1,nvwlf) = ia
                     ivwlf(2,nvwlf) = ic
c
c                 geometric combination rule
c
                     avwlf(nvwlf) = SQRT(avw_a * avw_c)
                     bvwlf(nvwlf) = SQRT(bve_a * bve_c)
                     if (mvw_a .eq. mvw_c) then
                        mvwlf(nvwlf) = mvw_a
                     else
                        write (iout,95) ita,itc
  95                    format (/,' KLFMM  --  Cannot process',
     &                       ' van der Waals-like L-L repulsion',/,
     &                       ' parameters having different exponent.',
     &                       /,' Check parameters for classes ',i6,
     &                       ' and ', i6)
                        call fatal
                     endif
                     donevw = .true.
                  endif
                  if (.not. donell .and. .not. donevw) then
                     write (iout,100) ita,itc
  100                format (/,' KLFMM  --  Unable to find a',
     &                    ' combination of parameters for',
     &                    ' ligand-ligand',/,
     &                    ' repulsion term involving ligand of class ',
     &                    i4,' and ',i4,';',/,
     &                    ' Check input and force field.')
                     call fatal
                  endif
               enddo
            enddo
         enddo
      endif
c
c     assign parameters for angular overlap model: calculation of lfse
c
      if (use_lfmmlfs) then
         do i = 1, nlfmm
            im = metlfmm(i)
            itm = class(im)
            size = 4
            call numeral (itm,pm,size)
            nlfse = nlfse + 1
            imlfse(nlfse) = im
            nelelfse(nlfse) = elelfmm(i)
            smlfse(nlfse) = spinlfmm(i)
            nl = 0
            if (n12(im) .eq. 0) then
               write(iout,105)
  105          format (/,' KLFMM  --  Metal atom is not connected to',
     &                 ' any ligand; Add metal-ligand bonds')
               call fatal
            endif
            do j = 1, n12(im)
               ia = i12(j,im)
               ita = class(ia)
               call numeral (ita,pa,size)
               if (ita .le. itm) then
                  pt = pa//pm
               else
                  pt = pm//pa
               end if
               foundsig = .false.
               foundpix = .false.
               foundpiy = .false.
               foundxds = .false.
               do k =1, maxnaom
                  if (aomsigid(k) .eq. pt) then
                     nl = nl +1
                     illfse(nlfse,nl) = ia
                     esig(nlfse,nl,1) = aomsig(k,1)
                     esig(nlfse,nl,2) = aomsig(k,2)
                     esig(nlfse,nl,3) = aomsig(k,3)
                     esig(nlfse,nl,4) = aomsig(k,4)
                     esig(nlfse,nl,5) = aomsig(k,5)
                     esig(nlfse,nl,6) = aomsig(k,6)
                     esig(nlfse,nl,7) = aomsig(k,7)
                     foundsig = .true.
                     exit
                  end if
               enddo
               if (.not. foundsig) then
                  write (iout,110) 'esig',itm,ita
  110             format (/,' KLFMM  --  Parameres not found for',
     &                    a4,' term to AOM (classes ',i5,' ',i5,').',
     &                    /,' Check input and force field file.')
                  call fatal
               endif
               do k =1, maxnaom
                  if (aompixid(k) .eq. pt) then
                     epix(nlfse,nl,1) = aompix(k,1)
                     epix(nlfse,nl,2) = aompix(k,2)
                     epix(nlfse,nl,3) = aompix(k,3)
                     epix(nlfse,nl,4) = aompix(k,4)
                     epix(nlfse,nl,5) = aompix(k,5)
                     epix(nlfse,nl,6) = aompix(k,6)
                     epix(nlfse,nl,7) = aompix(k,7)
                     foundpix = .true.
                     exit
                  end if
               enddo
               if (.not. foundpix) then
                  write (iout,110) 'epix',itm,ita
                  call fatal
               endif
               do k =1, maxnaom
                  if (aompiyid(k) .eq. pt) then
                     epiy(nlfse,nl,1) = aompiy(k,1)
                     epiy(nlfse,nl,2) = aompiy(k,2)
                     epiy(nlfse,nl,3) = aompiy(k,3)
                     epiy(nlfse,nl,4) = aompiy(k,4)
                     epiy(nlfse,nl,5) = aompiy(k,5)
                     epiy(nlfse,nl,6) = aompiy(k,6)
                     epiy(nlfse,nl,7) = aompiy(k,7)
                     foundpiy = .true.
                     exit
                  end if
               enddo
               if (.not. foundpiy) then
                  write (iout,110) 'epiy',itm,ita
                  call fatal
               endif 
               do k = 1, maxnaom
                  if (aomxdsid(k) .eq. pt) then
                     exds(nlfse,nl,1) = aomxds(k,1)
                     exds(nlfse,nl,2) = aomxds(k,2)
                     exds(nlfse,nl,3) = aomxds(k,3)
                     exds(nlfse,nl,4) = aomxds(k,4)
                     exds(nlfse,nl,5) = aomxds(k,5)
                     exds(nlfse,nl,6) = aomxds(k,6)
                     exds(nlfse,nl,7) = aomxds(k,7)
                     foundxds = .true.
                     exit
                  end if
               enddo
               if (.not. foundxds) then
                  write (iout,110) 'exds',itm,ita
                  call fatal
               endif
c
c           store unique ligand types
c
               do k = 1, maxlig
                  if (ligunqtyps(nlfse,k,1) .eq. ita) then
                     ligtyp(nlfse,nl) = k
                     exit
                  endif
                  if (ligunqtyps(nlfse,k,1) .eq. -1) then
                     ligunqtyps(nlfse,k,1) = ita
                     ligunqtyps(nlfse,k,2) = nl
                     ligtyp(nlfse,nl) = k
                     exit
                  endif
               enddo
            enddo
            nligs(nlfse) = nl
c
c        preparing reordered list of atoms in this aom region
c
            call reorderaom (nlfse)
         enddo
      endif
c
c     assign parameters for electron pairing energy model
c
      if (use_lfmmelp) then
         do i = 1, nlfmm
c
c        Set pairing energy only for low and intermediate spin states
c
            wtmp = 0.0d0
            if (spinlfmm(i) .eq. 0) then
               wtmp = 1.0d0
            else if (spinlfmm(i) .eq. 2) then
c
c TODO: set wtmp according to the number of d-electrons and adapt LFSE routine
c
  115          format (/,' KLFMM  --  Treatment of intermediate spin'
     &                 ' models not yet implemented. ')
               write (iout,115)
               call fatal
            else 
               cycle
            endif
            im = metlfmm(i)
            itm = class(im)
            size = 4
            call numeral (itm,pm,size)
            do j = 1, n12(im)
               ia = i12(j,im)
               ita = class(ia)
               call numeral (ita,pa,size)
               if (ita .le. itm) then
                  pt = pa//pm
               else
                  pt = pm//pa
               end if
               foundpair = .false.
               do k = 1, maxnaom
                  if (epairid(k) .eq. pt) then
                     npair = npair + 1
                     ipair(1,npair) = im
                     ipair(2,npair) = ia
                     wpair(npair) = wtmp
                     ppair(npair,1) = epairpar(k,1)
                     ppair(npair,2) = epairpar(k,2)
                     ppair(npair,3) = epairpar(k,3)
                     ppair(npair,4) = epairpar(k,4)
                     ppair(npair,5) = epairpar(k,5)
                     ppair(npair,6) = epairpar(k,6)
                     ppair(npair,7) = epairpar(k,7)
                     foundpair = .true.
                  end if
               end do
               if (.not. foundpair) then
                  write (iout,120) itm,ita
  120             format (/,' KLFMM  --  Unable to find LFMM parameters'
     &                    ' for electron pairing energy',/,
     &                    ' for classes ',i5,' and ',i5,');',/,
     &                    ' Check input and force field file.')
                  call fatal
               endif
            enddo
         enddo    
      endif
      return
      end
c
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine reorderaom  --  prepare reordered list of atoms for AOM  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "reorderaom" prepare a reordered list of atoms to satisfy the
c     needs of the LFSE routine for LFMM calculations. For a metal center
c     it defines a vector of atom indices according to this order:
c     1st block: METAL; 
c     2nd block: COORDINATING ATOMS;
c     3rd block: SUBSTIDIARY LIGAND atoms connected with 2nd block.
c     All the remaining atoms are ignored. The reordering list of 
c     atoms can be accesses per each LFMM center via the vectors
c     'aom2tnk' and 'tnk2aom' (see lfmm.i).
c
c
      subroutine reorderaom (ns)
      implicit  none
      include 'sizes.i'
      include 'couple.i'
      include 'lfmm.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,ii,ns
      integer il,ic,nssloc
      integer idmetal,nnotmet
      integer tnk2aom
      real*8 angle
c
c
c     initialization: nothing in AOM region
c
      totaom(ns) = 0
      do i = 1, maxcrd
         aom2tnk(ns,i) = 0
      enddo
      do i = 1, maxlig
         iss1(ns,i) = 0
         iss2(ns,i) = 0
         nss1(ns,i) = 0
         nss2(ns,i) = 0
      enddo
c
c     1st block: metal atom
c
      idmetal = imlfse(ns)
      aom2tnk(ns,1) = idmetal
      totaom(ns) = 1
      ii = 2  ! this is the index of the entry where to put the next atom
c
c     2nd block: coordinating atoms
c
      do i = 1, nligs(ns)
         aom2tnk(ns,ii) = illfse(ns,i)
         totaom(ns) = totaom(ns) + 1
         ii = ii + 1
      enddo
c
c     3rd block: the two atoms bounded to coordinating atoms (apart from metal)
c
      do i = 1, nligs(ns)
         il = illfse(ns,i)
         nnotmet = n12(il) - 1
         nssloc = 0
         if (nnotmet .eq. 0) then
            cycle
         elseif (nnotmet .ge. 3) then
            cycle
         elseif (nnotmet .eq. 1) then
            do j = 1, n12(il)
               ic = i12(j,il)
               if (ic .eq. idmetal) then
                  cycle
               else
                  call getangle (idmetal,il,ic,angle)
                  if (abs(angle - 180.0d0) .lt. 0.001 ) cycle
                  nssloc = nssloc + 1
                  nss1(ns,i) = 1
                  iss1(ns,i) = ic
               endif
            enddo
         elseif (nnotmet .eq. 2) then
            do j = 1, n12(il)
               ic = i12(j,il)
               if (ic .eq. idmetal) then
                  cycle
               else
                  nssloc = nssloc + 1
                  if (nssloc .eq. 1) then
                     nss1(ns,i) = 1
                     iss1(ns,i) = ic
                  else 
                     nss2(ns,i) = 1
                     iss2(ns,i) = ic
                  endif
               endif
            enddo
         endif
      enddo
      do i = 1, nligs(ns)
         if (nss1(ns,i) .eq. 1) then
            aom2tnk(ns,ii) = iss1(ns,i)
            totaom(ns) = totaom(ns) + 1
            ii = ii + 1
         endif
      enddo
      do i = 1, nligs(ns)
         if (nss2(ns,i) .eq. 1) then
            aom2tnk(ns,ii) = iss2(ns,i)
            totaom(ns) = totaom(ns) + 1
            ii = ii + 1
         endif
      enddo
      do i = 1, nligs(ns)
         if ((nss2(ns,i) .gt. 1) .or. (nss2(ns,i) .gt. 1)) then
            write(iout,20)
  20        format (/,' REORDERAOM  --  Unexpected ligand requiring'
     &              ' more than two second shell atoms. This case'
     &              ' need to be nefined in the code.')
            call fatal
         endif
      enddo
c
c     print information for debug    
c
      if (debug) then
        write(iout,30) ns,imlfse(ns)
  30    format(/,' REORDERAOM  --  AOM atoms list for LFMM center ',i3,
     &          ' (metal atom ',i6,')')
        do i = 1, totaom(ns)
          write(iout,40) aom2tnk(ns,i),tnk2aom (ns,aom2tnk(ns,i))
  40      format(5x,'atom ',i6,' has AOM index ',i3)
        enddo
        write(iout,50) ns,imlfse(ns)
  50    format(/,' REORDERAOM  -- ligands for LFMM center ',i3,
     &          ' (metal atom ',i6,')')
        do i = 1, nligs(ns)
          write(iout,60) i,illfse(ns,i)
  60      format (4x,'Ligand ',i2,' is atom ',i6)
          if (nss1(ns,i) .ne. 0) then
            write(iout,70) 1, iss1(ns,i)
  70        format(8x,'Subsidiary atom ',i1,': ',i6)
          endif
          if (nss2(ns,i) .ne. 0) then
            write(iout,70) 2, iss2(ns,i)
          endif
        enddo
      endif
      end subroutine
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function tnk2aom  --  get the atom index in AOM  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "tnk2aom" return the atom index in the reordered list of atoms:
c     n = tnk2aom(i,m) where m-th atom in TINKER list is the n-th
c     in the AOM region of LFMM center i.
c
c
      function tnk2aom (i,m)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'lfmm.i'
      integer i,j,m
      integer tnk2aom

      do j=1,totaom(i)
         if (aom2tnk(i,j) .eq. m) then
           tnk2aom=j
           return
         endif
      enddo
      write (iout,10) m,i
  10  format (/'  TNK2AOM  --  Atom ',i4,
     &        ' is not in aom region of center ',i4)
      call fatal
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine getangle  --  returns the angle in deg  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "getangle" calculates the value of the angle described by the
c     three atoms given as argument.
c
c
      subroutine getangle (ia,ib,ic,angle)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'math.i'
      integer ia,ib,ic 
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rcb2
      real*8 angle,dot
      real*8 cosine

      xia = x(ia)
      yia = y(ia)
      zia = z(ia)
      xib = x(ib)
      yib = y(ib)
      zib = z(ib)
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xab = xia - xib
      yab = yia - yib
      zab = zia - zib
      xcb = xic - xib
      ycb = yic - yib
      zcb = zic - zib
      rab2 = xab*xab + yab*yab + zab*zab
      rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
      dot = xab*xcb + yab*ycb + zab*zcb
      cosine = dot / sqrt(rab2*rcb2)
      cosine = min(1.0d0,max(-1.0d0,cosine))
      angle = radian * acos(cosine)

      return
      end subroutine
