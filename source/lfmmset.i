c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  lfmmset.i  --  general settings for LFMM calculation   ##
c     ##                                                         ##
c     #############################################################
c
c 
c     reproduceold logical flag behave as pre-2014 LFMM
c
c     shajorlfse   flag to use Shaffer and Jorgensen d/s-mixing
c     barysjlfse   barycenter of Shaffer and Jorgensen d/s-mixing
c     anlderlfse   flag to use analytic derivative  
c     stpzlfse     step size for finite difference
c     displfse     displacement for finite difference
c     doderlfse    flag controlling calculation of first derivatives
c     scnderlfse   flag controlling calculation of second derivatives
c     sdlfsedone   logical flag reporting status of LFSE second derivatives
c     dosdfdout    logical flag use finite diff lfse gradient to get Hessian
c     updtlfshess  logical flag uptate lfse contribution to Hessian
c     stbenglfse   flag controlling stabilization energy formula MOSE/LFSE
c     debuglfse    flag controlling debug printing of LFSE executable
c 
c     nlfmm        total number of lfmm centers
c
c     metlfmm      index of metal atom treated by LFMM per each lfmm center
c     elelfmm      number of electrons per each lfmm center
c     spinlfmm     spinstate flag (0:low; 1:high; 2:intermediate) per each lfmm center
c
c
      integer shajorlfse,anlderlfse,doderlfse
      integer debuglfse,scnderlfse,stbenglfse
      integer nlfmm,barysjlfse,displfse
      integer metlfmm,elelfmm,spinlfmm
      real*8 stpzlfse
      logical reproduceold,updtlfshess
      logical sdlfsedone,dosdfdout
      common /lfmmset/ stpzlfse,nlfmm,shajorlfse,doderlfse,
     &                 anlderlfse,debuglfse,scnderlfse,
     &                 stbenglfse,barysjlfse,displfse,
     &                 metlfmm(maxlfmm),elelfmm(maxlfmm),
     &                 spinlfmm(maxlfmm),
     &                 reproduceold,dosdfdout,
     &                 sdlfsedone(maxlfmm),updtlfshess
