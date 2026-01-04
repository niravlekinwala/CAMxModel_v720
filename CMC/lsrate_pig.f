      subroutine lsrate_pig(neq,t,y,rate,nr,r)
      implicit none
c
c----CAMx v7.20 220430
c
c     LSRATE_pig computes double precision species rates
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6
c
c     Routines Called:
c        LSRXN_pig
c
c     Called by:
c        LSODE
c
      include "camx.prm"
      include "chmdat.inc"
      include "lsbox.inc"
      include "ddmchm.inc"
c
      integer neq, nr, l
      double precision t
      double precision y(*)
      double precision rate(neq)
      double precision r(MXREACT)
      double precision loss(MXSPEC+1)
      double precision gain(MXSPEC+1)
c
c --- Entry point
c
c --- Get reaction rates, r
c
      call lsrxn_pig(dH2O,dM,dO2,dCH4,dH2,y,dbrk,r)
c
      do l = 1,neq
        Loss(l) = 0.0d0
        Gain(l) = 0.0d0
      enddo
c
c --- Calculate the species rates
c
c
c   NO2    NO     O    O3  HONO   O1D    OH  HNO3   NO3  N2O5
c   SO2  SULF   HO2    CO  FORM
c
        Loss(l_NO2 )= +         r(  1)+         r(  5)+         r( 10)
     &                +         r( 11)+         r( 15)+         r( 16)
c
        Gain(l_NO2 )= +         r(  3)+( 2.000)*r(  4)+         r( 12)
     &                +( 2.000)*r( 14)+         r( 15)+         r( 17)
     &                +         r( 23)
c
        Loss(l_NO  )= +         r(  3)+( 2.000)*r(  4)+         r(  5)
     &                +         r( 14)+         r( 23)
c
        Gain(l_NO  )= +         r(  1)+         r(  9)+         r( 13)
     &                +         r( 15)
c
        Loss(l_O   )= +         r(  2)
c
        Gain(l_O   )= +         r(  1)+         r(  7)+         r( 12)
c
        Loss(l_O3  )= +         r(  3)+         r(  6)+         r( 11)
c
        Gain(l_O3  )= +         r(  2)
c
        Loss(l_HONO)= +         r(  9)
c
        Gain(l_HONO)= +( 2.000)*r(  5)
c
        Loss(l_O1D )= +         r(  7)+         r(  8)
c
        Gain(l_O1D )= +         r(  6)
c
        Loss(l_OH  )= +         r( 10)+         r( 19)+         r( 20)
c
        Gain(l_OH  )= +( 2.000)*r(  8)+         r(  9)+         r( 23)
c
        Loss(l_HNO3)= 0.0
c
        Gain(l_HNO3)= +         r( 10)+( 2.000)*r( 18)
c
        Loss(l_NO3 )= +         r( 12)+         r( 13)+         r( 14)
     &                +         r( 15)+         r( 16)
c
        Gain(l_NO3 )= +         r( 11)+         r( 17)
c
        Loss(l_N2O5)= +         r( 17)+         r( 18)
c
        Gain(l_N2O5)= +         r( 16)
c
        Loss(l_SO2 )= +         r( 19)
c
        Gain(l_SO2 )= 0.0
c
        Loss(l_SULF)= 0.0
c
        Gain(l_SULF)= +         r( 19)
c
        Loss(l_HO2 )= +         r( 23)
c
        Gain(l_HO2 )= +         r( 19)+         r( 20)+( 2.000)*r( 21)
c
        Loss(l_CO  )= +         r( 20)
c
        Gain(l_CO  )= +         r( 21)+         r( 22)
c
        Loss(l_FORM)= +         r( 21)+         r( 22)
c
        Gain(l_FORM)= 0.0
c
c
      do l = 1,neq
        rate(l) = gain(l) -loss(l)
      enddo
c
      return
      end
