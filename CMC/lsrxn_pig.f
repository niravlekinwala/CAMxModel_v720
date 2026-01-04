      subroutine lsrxn_pig(H2O,M,O2,CH4,H2,y,dbrk,r)
      implicit none
c
c----CAMx v7.20 220430
c
c     LSRXN_pig computes double precision fluxes for each reaction
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6
c
c     Routines Called:
c        none
c
c     Called by:
c        LSRATE_pig
c
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
      double precision H2O, M, O2, CH4, H2, N2
      double precision y(*)
      double precision dbrk(*)
      double precision r(*)
c
c --- Entry point
c
      N2  = M - O2
      r(  1) = dbrk(  1)*y(l_NO2)
      r(  2) = dbrk(  2)*y(l_O)*O2*M
      r(  3) = dbrk(  3)*y(l_O3)*y(l_NO)
      r(  4) = dbrk(  4)*y(l_NO)*y(l_NO)*O2
      r(  5) = dbrk(  5)*y(l_NO)*y(l_NO2)*H2O
      r(  6) = dbrk(  6)*y(l_O3)
      r(  7) = dbrk(  7)*y(l_O1D)*M
      r(  8) = dbrk(  8)*y(l_O1D)*H2O
      r(  9) = dbrk(  9)*y(l_HONO)
      r( 10) = dbrk( 10)*y(l_NO2)*y(l_OH)
      r( 11) = dbrk( 11)*y(l_NO2)*y(l_O3)
      r( 12) = dbrk( 12)*y(l_NO3)
      r( 13) = dbrk( 13)*y(l_NO3)
      r( 14) = dbrk( 14)*y(l_NO3)*y(l_NO)
      r( 15) = dbrk( 15)*y(l_NO3)*y(l_NO2)
      r( 16) = dbrk( 16)*y(l_NO3)*y(l_NO2)
      r( 17) = dbrk( 17)*y(l_N2O5)
      r( 18) = dbrk( 18)*y(l_N2O5)*H2O
      r( 19) = dbrk( 19)*y(l_SO2)*y(l_OH)
      r( 20) = dbrk( 20)*y(l_OH)*y(l_CO)
      r( 21) = dbrk( 21)*y(l_FORM)
      r( 22) = dbrk( 22)*y(l_FORM)
      r( 23) = dbrk( 23)*y(l_HO2)*y(l_NO)
c
      return
      end
