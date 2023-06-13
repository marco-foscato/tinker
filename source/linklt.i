c
c
c     #####################################################################
c     ##                                                                 ##
c     ##  linklt.i  -  storage of information shared by Tinker and LFSE  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c
c     nmetal  
c     xyz     vector of coordinates [x1,...,xn;y1,...,yn;z1,...,zn]
c     grad    first derivatives with respect to coordinated
c     
      integer nmetalllt
      character(len=8) txt

      common /linklt/ nmetalllt,txt

