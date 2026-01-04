      subroutine ddmrate6(ny,nr,r,gain,loss)
      implicit none
c
c----CAMx v7.20 220430
c
c     DDMRATE6 computes species production and loss rates
c     needed only for DDM rate constant sensitivity
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        EBISOLV
c
c --- Argument definitions:
c        ny   - dimension of gain and loss
c        nr   - dimension of r
c        r    - reaction rates (hr-1)
c        gain - species production (ppm/hr)
c        loss - species destruction (ppm/hr)
c
c --- Includes:
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
c --- Arguments:
      integer ny, nr
      real    loss(ny+1), gain(ny+1)
      real    r(nr)
c
c --- Entry Point:
c
c
c-----Calculate the species rates
c
c
c    NO   NO2    O3     O   O1D    OH   HO2  HONO   PNA   PAN
c  C2O3   NO3  N2O5
c
        Loss(lNO   )= +         r(  3)+         r(  6)+         r( 16)
     &                +( 2.000)*r( 22)+         r( 23)+         r( 24)
     &                +         r( 30)+         r( 54)+         r( 55)
     &                +         r( 68)+         r( 81)+         r( 88)
     &                +         r(103)+         r(132)
c
        Gain(lNO   )= +         r(  1)+         r(  4)+         r( 15)
     &                +         r( 17)+         r( 25)+         r( 27)
     &                +( 0.200)*r(148)
c
        Loss(lNO2  )= +         r(  1)+         r(  4)+         r(  5)
     &                +         r(  7)+         r( 17)+         r( 18)
     &                +         r( 23)+         r( 28)+         r( 31)
     &                +         r( 89)+         r(104)+         r(118)
     &                +         r(136)+         r(148)
c
        Gain(lNO2  )= +         r(  3)+         r(  6)+         r( 14)
     &                +( 2.000)*r( 16)+         r( 17)+         r( 21)
     &                +( 2.000)*r( 22)+         r( 26)+         r( 27)
     &                +         r( 30)+         r( 32)+         r( 33)
     &                +         r( 46)+         r( 47)+         r( 49)
     &                +( 2.000)*r( 50)+( 0.610)*r( 51)+         r( 52)
     &                +         r( 53)+         r( 54)+         r( 62)
     &                +         r( 68)+         r( 81)+         r( 88)
     &                +         r( 90)+         r( 91)+         r(103)
        Gain(lNO2  ) = Gain(lNO2  )
     &                +         r(105)+         r(106)+         r(107)
     &                +         r(122)+         r(126)+         r(130)
     &                +( 0.900)*r(132)+( 0.200)*r(147)+( 0.470)*r(156)
c
        Loss(lO3   )= +         r(  3)+         r(  7)+         r(  8)
     &                +         r(  9)+         r( 12)+         r( 13)
     &                +         r( 49)+         r(121)+         r(125)
     &                +         r(129)+         r(140)+         r(146)
     &                +         r(150)+         r(155)
c
        Gain(lO3   )= +         r(  2)+( 0.200)*r( 92)+( 0.200)*r(108)
c
        Loss(lO    )= +         r(  2)+         r(  4)+         r(  5)
     &                +         r(  6)+         r( 40)+         r( 44)
     &                +         r( 45)+         r( 46)+         r( 77)
     &                +         r( 84)+         r( 99)+         r(119)
     &                +         r(123)+         r(127)+         r(144)
     &                +         r(153)
c
        Gain(lO    )= +         r(  1)+         r(  8)+         r( 10)
     &                +         r( 14)+         r( 41)+( 0.500)*r(129)
c
        Loss(lO1D  )= +         r( 10)+         r( 11)+         r( 38)
c
        Gain(lO1D  )= +         r(  9)
c
        Loss(lOH   )= +         r( 12)+         r( 24)+         r( 26)
     &                +         r( 28)+         r( 29)+         r( 33)
     &                +         r( 37)+         r( 39)+         r( 40)
     &                +( 2.000)*r( 41)+( 2.000)*r( 42)+         r( 43)
     &                +         r( 47)+         r( 61)+         r( 63)
     &                +         r( 64)+         r( 66)+         r( 67)
     &                +         r( 71)+         r( 73)+         r( 74)
     &                +         r( 83)+         r( 85)+         r( 96)
     &                +         r( 98)+         r(100)+         r(107)
        Loss(lOH   ) = Loss(lOH   )
     &                +         r(113)+         r(114)+         r(115)
     &                +         r(120)+         r(124)+         r(128)
     &                +         r(131)+         r(134)+         r(139)
     &                +         r(141)+         r(142)+         r(145)
     &                +         r(149)+         r(154)
c
        Gain(lOH   )= +( 2.000)*r( 11)+         r( 13)+         r( 25)
     &                +         r( 30)+( 2.000)*r( 36)+         r( 38)
     &                +         r( 44)+         r( 45)+( 0.390)*r( 51)
     &                +         r( 52)+         r( 65)+         r( 72)
     &                +         r( 77)+         r( 84)+         r( 97)
     &                +         r( 99)+( 0.100)*r(119)+( 0.100)*r(121)
     &                +( 0.300)*r(123)+( 0.130)*r(125)+( 0.500)*r(129)
     &                +( 0.080)*r(140)+( 0.266)*r(146)+( 0.268)*r(150)
     &                +( 0.570)*r(155)
c
        Loss(lHO2  )= +         r( 13)+         r( 30)+         r( 31)
     &                +( 2.000)*r( 34)+( 2.000)*r( 35)+         r( 43)
     &                +         r( 44)+         r( 48)+         r( 56)
     &                +         r( 57)+         r( 69)+         r( 79)
     &                +         r( 82)+         r( 92)+         r(108)
     &                +         r(137)
c
        Gain(lHO2  )= +         r( 12)+         r( 32)+         r( 37)
     &                +         r( 38)+         r( 39)+         r( 40)
     &                +         r( 45)+         r( 47)+( 0.610)*r( 51)
     &                +         r( 61)+         r( 62)+         r( 63)
     &                +         r( 65)+         r( 66)+         r( 68)
     &                +( 0.740)*r( 70)+( 0.300)*r( 71)+         r( 72)
     &                +         r( 73)+         r( 74)+( 2.000)*r( 75)
     &                +         r( 77)+         r( 78)+         r( 80)
     &                +         r( 81)+         r( 83)+         r( 87)
        Gain(lHO2  ) = Gain(lHO2  )
     &                +( 0.900)*r( 93)+         r(102)+         r(103)
     &                +         r(109)+( 2.000)*r(111)+         r(112)
     &                +         r(113)+         r(114)+( 0.110)*r(115)
     &                +( 0.940)*r(116)+         r(117)+( 0.300)*r(119)
     &                +( 0.950)*r(120)+( 0.440)*r(121)+( 1.700)*r(123)
     &                +         r(124)+( 0.130)*r(125)+( 0.100)*r(127)
     &                +         r(128)+( 0.500)*r(129)+         r(130)
     &                +( 0.440)*r(131)+( 0.900)*r(132)+         r(133)
        Gain(lHO2  ) = Gain(lHO2  )
     &                +( 0.600)*r(134)+         r(138)+( 2.000)*r(139)
     &                +( 0.760)*r(140)+( 0.700)*r(141)+         r(143)
     &                +( 0.250)*r(144)+( 0.912)*r(145)+( 0.066)*r(146)
     &                +( 0.800)*r(147)+( 0.800)*r(148)+( 0.503)*r(149)
     &                +( 0.154)*r(150)+( 0.925)*r(151)+( 1.033)*r(152)
     &                +( 0.750)*r(154)+( 0.070)*r(155)+( 0.280)*r(156)
c
        Loss(lHONO )= +         r( 25)+         r( 26)+( 2.000)*r( 27)
c
        Gain(lHONO )= +( 2.000)*r( 23)+         r( 24)
c
        Loss(lPNA  )= +         r( 32)+         r( 33)+         r( 51)
c
        Gain(lPNA  )= +         r( 31)
c
        Loss(lPAN  )= +         r( 90)+         r( 91)
c
        Gain(lPAN  )= +         r( 89)
c
        Loss(lC2O3 )= +         r( 88)+         r( 89)+         r( 92)
     &                +         r( 93)+         r( 94)+( 2.000)*r( 95)
     &                +         r(112)
c
        Gain(lC2O3 )= +         r( 84)+         r( 85)+         r( 86)
     &                +         r( 90)+         r( 91)+         r( 96)
     &                +         r(138)+         r(139)+( 0.620)*r(140)
     &                +         r(142)+         r(143)+( 0.210)*r(149)
     &                +( 0.114)*r(150)+( 0.967)*r(152)
c
        Loss(lNO3  )= +         r( 14)+         r( 15)+         r( 16)
     &                +         r( 17)+         r( 18)+         r( 46)
     &                +         r( 47)+         r( 48)+         r( 49)
     &                +( 2.000)*r( 50)+         r( 78)+         r( 86)
     &                +         r(101)+         r(122)+         r(126)
     &                +         r(130)+         r(135)+         r(147)
     &                +         r(151)+         r(156)
c
        Gain(lNO3  )= +         r(  5)+         r(  7)+         r( 21)
     &                +         r( 29)+( 0.390)*r( 51)+         r( 53)
c
        Loss(lN2O5 )= +         r( 19)+         r( 20)+         r( 21)
     &                +         r( 53)
c
        Gain(lN2O5 )= +         r( 18)
c
c
      return
      end

