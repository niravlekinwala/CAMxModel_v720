      subroutine hr_xo3(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_XO3 solves halogen oxides by EBI iterations
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
c        y0  - initial Conc for this step (ppm)
c        yh  - current Conc iteration (ppm)
c        y1  - next Conc iteration (ppm)
c        H2O - water vapor Conc (ppm)
c        M   - total gas Conc (ppm)
c        O2  - oxygen Conc (ppm)
c        CH4 - methane Conc (ppm)
c        H2  - hydrogen Conc (ppm)
c        ny  - dimension of y0, y1 and yh
c        rk  - rate constants (ppm-n hr-1)
c        r   - reaction rates (hr-1)
c        nr  - dimension of rk and r
c        dt  - time step (hr)
c
c --- Includes:
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
c --- Arguments:
      integer ny, nr
      real y0(ny+1), y1(ny+1), yh(ny+1)
      real rk(nr), r(nr)
      real H2O, M, O2, CH4, H2, dt
c
c --- Local variables:
      real N2
      integer iter
      real TINY
      parameter (TINY = 1.0E-15)
      real gainCL, lossCL
      real gainCLO, lossCLO
      real gainBR, lossBR
      real gainBRO, lossBRO
      real gainI, lossI
      real gainIO, lossIO
      real gainOIO, lossOIO
      real r217
      real r218
      real r219
      real r220
      real r221
      real r222
      real r223
      real r224
      real r225
      real r226
      real r228
      real r229
      real r230
      real r231
      real r232
      real r233
      real r234
      real r235
      real r236
      real r238
      real r239
      real r240
      real r241
      real r242
      real r243
      real r244
      real r245
      real r246
      real r247
      real r248
      real r249
      real r250
      real r251
      real r252
      real r253
      real r254
      real r255
      real r257
      real r258
      real r259
      real r260
      real r261
      real r262
      real r263
      real r264
      real r265
      real r266
      real r267
      real r268
      real r269
      real r270
      real r271
      real r272
      real r273
      real r274
      real r275
      real r276
      real r277
      real r278
      real r279
      real r280
      real r281
      real r282
      real r283
      real r285
      real r287
      real r289
      real r290
      real r291
      real r292
      real r293
      real r294
      real r295
      real r296
      real r297
      real r298
      real r299
      real r300
      real r301
      real r302
c
c --- Entry Point
c
      N2 = M - O2
c
c --- First iteration uses reaction rates from r()
c
      r217 = r(217)
      r218 = r(218)
      r219 = r(219)
      r220 = r(220)
      r221 = r(221)
      r222 = r(222)
      r223 = r(223)
      r224 = r(224)
      r225 = r(225)
      r226 = r(226)
      r228 = r(228)
      r229 = r(229)
      r230 = r(230)
      r231 = r(231)
      r232 = r(232)
      r233 = r(233)
      r234 = r(234)
      r235 = r(235)
      r236 = r(236)
      r238 = r(238)
      r239 = r(239)
      r240 = r(240)
      r241 = r(241)
      r242 = r(242)
      r243 = r(243)
      r244 = r(244)
      r245 = r(245)
      r246 = r(246)
      r247 = r(247)
      r248 = r(248)
      r249 = r(249)
      r250 = r(250)
      r251 = r(251)
      r252 = r(252)
      r253 = r(253)
      r254 = r(254)
      r255 = r(255)
      r257 = r(257)
      r258 = r(258)
      r259 = r(259)
      r260 = r(260)
      r261 = r(261)
      r262 = r(262)
      r263 = r(263)
      r264 = r(264)
      r265 = r(265)
      r266 = r(266)
      r267 = r(267)
      r268 = r(268)
      r269 = r(269)
      r270 = r(270)
      r271 = r(271)
      r272 = r(272)
      r273 = r(273)
      r274 = r(274)
      r275 = r(275)
      r276 = r(276)
      r277 = r(277)
      r278 = r(278)
      r279 = r(279)
      r280 = r(280)
      r281 = r(281)
      r282 = r(282)
      r283 = r(283)
      r285 = r(285)
      r287 = r(287)
      r289 = r(289)
      r290 = r(290)
      r291 = r(291)
      r292 = r(292)
      r293 = r(293)
      r294 = r(294)
      r295 = r(295)
      r296 = r(296)
      r297 = r(297)
      r298 = r(298)
      r299 = r(299)
      r300 = r(300)
      r301 = r(301)
      r302 = r(302)
c
c --- Gain and loss terms for the XO species
c
      gainCL = 0.0
     &       + ( 2.000)*r217
     &       +          r218
     &       + ( 1.400)*r220
     &       +          r221
     &       +          r226
     &       +          r228
     &       +          r229
     &       +          r230
     &       +          r231
     &       +          r238
     &       +          r289
     &       +          r290
     &       +          r295
c
      lossCL = 0.0
     &       +          r219
     &       +          r232
     &       +          r233
     &       +          r234
     &       +          r235
     &       +          r236
c
      gainCLO = 0.0
     &       +          r219
     &       +          r224
     &       +          r225
c
      lossCLO = 0.0
     &       + ( 2.000)*r220
     &       +          r221
     &       +          r222
     &       +          r223
     &       +          r231
     &       +          r289
     &       +          r290
c
      gainBR = 0.0
     &       + ( 2.000)*r239
     &       +          r240
     &       +          r241
     &       +          r242
     &       +          r247
     &       +          r249
     &       + ( 2.000)*r250
     &       +          r252
     &       +          r254
     &       +          r255
     &       +          r257
     &       +          r258
     &       + ( 0.250)*r259
     &       +          r289
     &       +          r291
     &       +          r294
     &       + ( 3.000)*r296
     &       + ( 3.000)*r297
     &       + ( 2.000)*r298
     &       +          r299
     &       +          r300
     &       +          r301
c
      lossBR = 0.0
     &       +          r243
     &       +          r244
     &       +          r245
     &       +          r246
     &       +          r260
     &       +          r261
     &       +          r262
     &       +          r263
c
      gainBRO = 0.0
     &       +          r243
     &       +          r246
c
      lossBRO = 0.0
     &       +          r247
     &       +          r248
     &       +          r249
     &       + ( 2.000)*r250
     &       + ( 2.000)*r251
     &       +          r252
     &       +          r253
     &       +          r259
     &       +          r289
     &       +          r291
c
      gainI = 0.0
     &       + ( 2.000)*r264
     &       +          r265
     &       +          r266
     &       +          r267
     &       +          r268
     &       +          r272
     &       + ( 0.400)*r273
     &       +          r275
     &       +          r278
     &       +          r283
     &       +          r285
     &       +          r287
     &       +          r290
     &       +          r291
     &       +          r292
     &       + ( 2.000)*r293
     &       +          r294
     &       +          r295
     &       +          r302
c
      lossI = 0.0
     &       +          r269
     &       +          r270
     &       +          r271
c
      gainIO = 0.0
     &       +          r269
     &       +          r277
     &       +          r282
c
      lossIO = 0.0
     &       +          r272
     &       + ( 2.000)*r273
     &       +          r274
     &       +          r275
     &       +          r276
     &       +          r280
     &       +          r290
     &       +          r291
c
      gainOIO = 0.0
     &       + ( 0.400)*r273
     &       +          r283
c
      lossOIO = 0.0
     &       +          r278
     &       +          r279
     &       +          r280
     &       + ( 2.000)*r281
     &       +          r282
c
c
c --- EBI solution for the XO species
c
      y1(lCL) = MIN(1.0, MAX(TINY,(y0(lCL) + gainCL*dt)
     &                         / (1.0 + lossCL*dt/yh(lCL))))
c
      y1(lCLO) = MIN(1.0, MAX(TINY,(y0(lCLO) + gainCLO*dt)
     &                         / (1.0 + lossCLO*dt/yh(lCLO))))
c
      y1(lBR) = MIN(1.0, MAX(TINY,(y0(lBR) + gainBR*dt)
     &                         / (1.0 + lossBR*dt/yh(lBR))))
c
      y1(lBRO) = MIN(1.0, MAX(TINY,(y0(lBRO) + gainBRO*dt)
     &                         / (1.0 + lossBRO*dt/yh(lBRO))))
c
      y1(lI) = MIN(1.0, MAX(TINY,(y0(lI) + gainI*dt)
     &                         / (1.0 + lossI*dt/yh(lI))))
c
      y1(lIO) = MIN(1.0, MAX(TINY,(y0(lIO) + gainIO*dt)
     &                         / (1.0 + lossIO*dt/yh(lIO))))
c
      y1(lOIO) = MIN(1.0, MAX(TINY,(y0(lOIO) + gainOIO*dt)
     &                         / (1.0 + lossOIO*dt/yh(lOIO))))
c
c --- Perform iteratations > 1
c     convergence will be tested in EBISOLV
c
      do iter = 2,3
c
c --- Update reaction rates that change with each iteration
c
      r219 = rk(219)*y1(lCL)*yh(lO3)
      r220 = rk(220)*y1(lCLO)*y1(lCLO)
      r221 = rk(221)*y1(lCLO)*yh(lNO)
      r222 = rk(222)*y1(lCLO)*yh(lHO2)
      r223 = rk(223)*y1(lCLO)*yh(lNO2)
      r231 = rk(231)*y1(lCLO)*yh(lMEO2)
      r232 = rk(232)*y1(lCL)*CH4
      r233 = rk(233)*y1(lCL)*yh(lPAR)
      r234 = rk(234)*y1(lCL)*yh(lETHA)
      r235 = rk(235)*y1(lCL)*yh(lPRPA)
      r236 = rk(236)*y1(lCL)*yh(lISOP)
      r243 = rk(243)*y1(lBR)*yh(lO3)
      r244 = rk(244)*y1(lBR)*yh(lHO2)
      r245 = rk(245)*y1(lBR)*yh(lNO2)
      r246 = rk(246)*y1(lBR)*yh(lNO3)
      r247 = rk(247)*y1(lBRO)
      r248 = rk(248)*y1(lBRO)*yh(lHO2)
      r249 = rk(249)*y1(lBRO)*yh(lOH)
      r250 = rk(250)*y1(lBRO)*y1(lBRO)
      r251 = rk(251)*y1(lBRO)*y1(lBRO)
      r252 = rk(252)*y1(lBRO)*yh(lNO)
      r253 = rk(253)*y1(lBRO)*yh(lNO2)
      r259 = rk(259)*y1(lBRO)*yh(lMEO2)
      r260 = rk(260)*y1(lBR)*yh(lFORM)
      r261 = rk(261)*y1(lBR)*yh(lALD2)
      r262 = rk(262)*y1(lBR)*yh(lOLE)
      r263 = rk(263)*y1(lBR)*yh(lISOP)
      r269 = rk(269)*y1(lI)*yh(lO3)
      r270 = rk(270)*y1(lI)*yh(lHO2)
      r271 = rk(271)*y1(lI)*yh(lNO2)
      r272 = rk(272)*y1(lIO)
      r273 = rk(273)*y1(lIO)*y1(lIO)
      r274 = rk(274)*y1(lIO)*yh(lHO2)
      r275 = rk(275)*y1(lIO)*yh(lNO)
      r276 = rk(276)*y1(lIO)*yh(lNO2)
      r278 = rk(278)*y1(lOIO)
      r279 = rk(279)*y1(lOIO)*yh(lOH)
      r280 = rk(280)*y1(lOIO)*y1(lIO)
      r281 = rk(281)*y1(lOIO)*y1(lOIO)
      r282 = rk(282)*y1(lOIO)*yh(lNO)
      r289 = rk(289)*y1(lCLO)*y1(lBRO)
      r290 = rk(290)*y1(lCLO)*y1(lIO)
      r291 = rk(291)*y1(lBRO)*y1(lIO)
c
c --- Gain and loss terms for the XO species
c
      gainCL = 0.0
     &       + ( 2.000)*r217
     &       +          r218
     &       + ( 1.400)*r220
     &       +          r221
     &       +          r226
     &       +          r228
     &       +          r229
     &       +          r230
     &       +          r231
     &       +          r238
     &       +          r289
     &       +          r290
     &       +          r295
c
      lossCL = 0.0
     &       +          r219
     &       +          r232
     &       +          r233
     &       +          r234
     &       +          r235
     &       +          r236
c
      gainCLO = 0.0
     &       +          r219
     &       +          r224
     &       +          r225
c
      lossCLO = 0.0
     &       + ( 2.000)*r220
     &       +          r221
     &       +          r222
     &       +          r223
     &       +          r231
     &       +          r289
     &       +          r290
c
      gainBR = 0.0
     &       + ( 2.000)*r239
     &       +          r240
     &       +          r241
     &       +          r242
     &       +          r247
     &       +          r249
     &       + ( 2.000)*r250
     &       +          r252
     &       +          r254
     &       +          r255
     &       +          r257
     &       +          r258
     &       + ( 0.250)*r259
     &       +          r289
     &       +          r291
     &       +          r294
     &       + ( 3.000)*r296
     &       + ( 3.000)*r297
     &       + ( 2.000)*r298
     &       +          r299
     &       +          r300
     &       +          r301
c
      lossBR = 0.0
     &       +          r243
     &       +          r244
     &       +          r245
     &       +          r246
     &       +          r260
     &       +          r261
     &       +          r262
     &       +          r263
c
      gainBRO = 0.0
     &       +          r243
     &       +          r246
c
      lossBRO = 0.0
     &       +          r247
     &       +          r248
     &       +          r249
     &       + ( 2.000)*r250
     &       + ( 2.000)*r251
     &       +          r252
     &       +          r253
     &       +          r259
     &       +          r289
     &       +          r291
c
      gainI = 0.0
     &       + ( 2.000)*r264
     &       +          r265
     &       +          r266
     &       +          r267
     &       +          r268
     &       +          r272
     &       + ( 0.400)*r273
     &       +          r275
     &       +          r278
     &       +          r283
     &       +          r285
     &       +          r287
     &       +          r290
     &       +          r291
     &       +          r292
     &       + ( 2.000)*r293
     &       +          r294
     &       +          r295
     &       +          r302
c
      lossI = 0.0
     &       +          r269
     &       +          r270
     &       +          r271
c
      gainIO = 0.0
     &       +          r269
     &       +          r277
     &       +          r282
c
      lossIO = 0.0
     &       +          r272
     &       + ( 2.000)*r273
     &       +          r274
     &       +          r275
     &       +          r276
     &       +          r280
     &       +          r290
     &       +          r291
c
      gainOIO = 0.0
     &       + ( 0.400)*r273
     &       +          r283
c
      lossOIO = 0.0
     &       +          r278
     &       +          r279
     &       +          r280
     &       + ( 2.000)*r281
     &       +          r282
c
c --- EBI solution for the XO species
c
      y1(lCL) = MIN(1.0, MAX(TINY,(y0(lCL) + gainCL*dt)
     &                         / (1.0 + lossCL*dt/y1(lCL))))
c
      y1(lCLO) = MIN(1.0, MAX(TINY,(y0(lCLO) + gainCLO*dt)
     &                         / (1.0 + lossCLO*dt/y1(lCLO))))
c
      y1(lBR) = MIN(1.0, MAX(TINY,(y0(lBR) + gainBR*dt)
     &                         / (1.0 + lossBR*dt/y1(lBR))))
c
      y1(lBRO) = MIN(1.0, MAX(TINY,(y0(lBRO) + gainBRO*dt)
     &                         / (1.0 + lossBRO*dt/y1(lBRO))))
c
      y1(lI) = MIN(1.0, MAX(TINY,(y0(lI) + gainI*dt)
     &                         / (1.0 + lossI*dt/y1(lI))))
c
      y1(lIO) = MIN(1.0, MAX(TINY,(y0(lIO) + gainIO*dt)
     &                         / (1.0 + lossIO*dt/y1(lIO))))
c
      y1(lOIO) = MIN(1.0, MAX(TINY,(y0(lOIO) + gainOIO*dt)
     &                         / (1.0 + lossOIO*dt/y1(lOIO))))
c
      enddo
c
      return
      end

