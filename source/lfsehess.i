c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2014 by Marco Foscato & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  lfsehess.i  --  Cartesian Hessian elements from LFSE  ##
c     ##                                                        ##
c     ############################################################
c
c
c     hellfse  vectors of Hessian contributions from for LFSE per each LFMM center
c     hidlfse  identification strings for Hessian contributions from for LFSE per each LFMM center
c
c
      real*8 hellfse
      character*16 hidlfse
      common /hlfse/ hellfse(maxlfmm,maxhlfse),
     &               hidlfse(maxlfmm,maxhlfse)
