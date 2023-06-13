c
c
c     #################################################################
c     ##                                                             ##
c     ##  ksnbnd.i -- parameters for special 12-10 non-bonded terms  ##
c     ##                                                             ##
c     #################################################################
c
c
c
c     maxpsnb  maximum number of parameters for special non-bonded terms
c   
c     kisnb     atom classes involved in special non-bonded terms
c     krsnd     radius parameter for special non-bonded terms
c     kepssnb   well depth parameter for special non-bonded terms
c     
c
      integer maxpsnb
      parameter (maxpsnb=20)
      integer kisnb
      real*8 krsnd,kepssnb
      common /ksnbnd/ kisnb(2,maxnsnb), 
     &                krsnd(maxnsnb),
     &                kepssnb(maxnsnb)
