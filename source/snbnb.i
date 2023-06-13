c
c
c     ###################################################################
c     ##                                                               ##
c     ##  snbnd.i  -- 12-10 non-bonded terms parameters for structure  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     nsnb     total number of special non-bonded terms 
c     isnb     atom numbers involved in each special non-bonded term
c     rsnd     radius parameter for each special non-bonded term
c     epssnb   well depth parameter for each special non-bonded term
c
c
      integer nsnb,isnb
      real*8 rsnd,epssnb
      common /snbnd/ isnb(2,maxnsnb), 
     &               rsnd(maxnsnb),
     &               epssnb(maxnsnb),
     &               nsnb
