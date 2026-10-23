\ Model amplgurobi
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x330572cf3dc04480
Minimize
  0 x[4] + 0 x[5] + 0 x[6] + 0 Obj_flatv1_ + 0 Obj_flatv3_ + 0 Obj_flatv4_
   + 0 Obj_flatv5_ + 0 Obj_flatv6_ + 0 Obj_flatv7_ + Obj_auxv
Subject To
 C1: x[1] + 40 x[2] = 176
 Obj__quadfn2_: - 400 x[3] - Obj_flatv2_ + [ x[3] ^2 ] = -40000
Bounds
 x[1] <= 200
 -10 <= x[2] <= 10
 -5 <= x[3] <= 12
 -15 <= x[4] <= 20
 -100 <= x[5] <= 1
 -5 <= x[6] <= 22
 Obj_flatv1_ = 0
 Obj_flatv2_ free
 Obj_flatv3_ = 0
 Obj_flatv4_ = 0
 Obj_flatv5_ = 0
 Obj_flatv6_ = 0
 Obj_flatv7_ = 0
 Obj_auxv free
General Constraints
\ Obj_auxv = +(-1 * exp(-2.20002 + (0.04 * x[4]) + (0.0004 * x[5]) + (0.04 
\ * x[6]) + (x[4] * x[4] * -0.002) + (x[5] * x[5] * -0.002) + (x[6] * x[6] 
\ * -0.0002) + ((1.4021754535147394e-05 + (-1.7272534013605442e-06 * x[1]) 
\ + (-1.0770975056689343e-05 * x[2]) + (x[1] * x[1] * 
\ 3.3216411564625849e-08) + (x[1] * x[2] * 1.0629251700680272e-06) + (x[2] 
\ * x[2] * 7.0861678004535146e-08)) * Obj_flatv2_ * -0.5) + 
\ (((-4.2517006802721087e-02 * x[4]) + (x[3] * x[4] * 
\ 2.1258503401360543e-04)) * ((-0.003125 * x[1]) + 
\ (-5.0333333333333334e-02 * x[2]) + (x[1] * x[1] * 0.00015625) + (x[1] * 
\ x[2] * 0.005) + (x[2] * x[2] * 3.3333333333333332e-04) + (Obj_flatv2_ * 
\ (6.5958333333333341e-02 + (-0.055 * x[2]) + (-0.003125 * x[3]) + 
\ (-6.6666666666666664e-04 * x[4]) + (x[2] * x[2] * 0.005) + (x[3] * x[3] 
\ * 0.00015625) + (x[4] * x[4] * 3.3333333333333332e-04)) * 
\ 2.1258503401360543e-04)) * -0.5)))
 Obj__linfn: Obj_auxv = NL : ( PLUS , -1 , -1 ) ( MULTIPLY , -1 , 0 ) 
  ( CONSTANT , -1 , 1 ) ( EXP , -1 , 1 ) ( PLUS , -1 , 3 ) 
  ( CONSTANT , -2.20002 , 4 ) ( MULTIPLY , -1 , 4 ) 
  ( CONSTANT , 0.04 , 6 ) ( VARIABLE , x[4] , 6 ) ( MULTIPLY , -1 , 4 ) 
  ( CONSTANT , 0.0004 , 9 ) ( VARIABLE , x[5] , 9 ) ( MULTIPLY , -1 , 4 ) 
  ( CONSTANT , 0.04 , 12 ) ( VARIABLE , x[6] , 12 ) ( MULTIPLY , -1 , 4 ) 
  ( VARIABLE , x[4] , 15 ) ( VARIABLE , x[4] , 15 ) 
  ( CONSTANT , -0.002 , 15 ) ( MULTIPLY , -1 , 4 ) 
  ( VARIABLE , x[5] , 19 ) ( VARIABLE , x[5] , 19 ) 
  ( CONSTANT , -0.002 , 19 ) ( MULTIPLY , -1 , 4 ) 
  ( VARIABLE , x[6] , 23 ) ( VARIABLE , x[6] , 23 ) 
  ( CONSTANT , -0.0002 , 23 ) ( MULTIPLY , -1 , 4 ) ( PLUS , -1 , 27 ) 
  ( CONSTANT , 1.40218e-05 , 28 ) ( MULTIPLY , -1 , 28 ) 
  ( CONSTANT , -1.72725e-06 , 30 ) ( VARIABLE , x[1] , 30 ) 
  ( MULTIPLY , -1 , 28 ) ( CONSTANT , -1.0771e-05 , 33 ) 
  ( VARIABLE , x[2] , 33 ) ( MULTIPLY , -1 , 28 ) ( VARIABLE , x[1] , 36 ) 
  ( VARIABLE , x[1] , 36 ) ( CONSTANT , 3.32164e-08 , 36 ) 
  ( MULTIPLY , -1 , 28 ) ( VARIABLE , x[1] , 40 ) ( VARIABLE , x[2] , 40 ) 
  ( CONSTANT , 1.06293e-06 , 40 ) ( MULTIPLY , -1 , 28 ) 
  ( VARIABLE , x[2] , 44 ) ( VARIABLE , x[2] , 44 ) 
  ( CONSTANT , 7.08617e-08 , 44 ) ( VARIABLE , Obj_flatv2_ , 27 ) 
  ( CONSTANT , -0.5 , 27 ) ( MULTIPLY , -1 , 4 ) ( PLUS , -1 , 50 ) 
  ( MULTIPLY , -1 , 51 ) ( CONSTANT , -0.0425170068027211 , 52 ) 
  ( VARIABLE , x[4] , 52 ) ( MULTIPLY , -1 , 51 ) ( VARIABLE , x[3] , 55 ) 
  ( VARIABLE , x[4] , 55 ) ( CONSTANT , 0.000212585 , 55 ) 
  ( PLUS , -1 , 50 ) ( MULTIPLY , -1 , 59 ) ( CONSTANT , -0.003125 , 60 ) 
  ( VARIABLE , x[1] , 60 ) ( MULTIPLY , -1 , 59 ) 
  ( CONSTANT , -0.0503333333333333 , 63 ) ( VARIABLE , x[2] , 63 ) 
  ( MULTIPLY , -1 , 59 ) ( VARIABLE , x[1] , 66 ) ( VARIABLE , x[1] , 66 ) 
  ( CONSTANT , 0.00015625 , 66 ) ( MULTIPLY , -1 , 59 ) 
  ( VARIABLE , x[1] , 70 ) ( VARIABLE , x[2] , 70 ) 
  ( CONSTANT , 0.005 , 70 ) ( MULTIPLY , -1 , 59 ) 
  ( VARIABLE , x[2] , 74 ) ( VARIABLE , x[2] , 74 ) 
  ( CONSTANT , 3.3333333333333332e-04 , 74 ) ( MULTIPLY , -1 , 59 ) 
  ( VARIABLE , Obj_flatv2_ , 78 ) ( PLUS , -1 , 78 ) 
  ( CONSTANT , 0.0659583333333333 , 80 ) ( MULTIPLY , -1 , 80 ) 
  ( CONSTANT , -0.055 , 82 ) ( VARIABLE , x[2] , 82 ) 
  ( MULTIPLY , -1 , 80 ) ( CONSTANT , -0.003125 , 85 ) 
  ( VARIABLE , x[3] , 85 ) ( MULTIPLY , -1 , 80 ) 
  ( CONSTANT , -6.6666666666666664e-04 , 88 ) ( VARIABLE , x[4] , 88 ) 
  ( MULTIPLY , -1 , 80 ) ( VARIABLE , x[2] , 91 ) ( VARIABLE , x[2] , 91 ) 
  ( CONSTANT , 0.005 , 91 ) ( MULTIPLY , -1 , 80 ) 
  ( VARIABLE , x[3] , 95 ) ( VARIABLE , x[3] , 95 ) 
  ( CONSTANT , 0.00015625 , 95 ) ( MULTIPLY , -1 , 80 ) 
  ( VARIABLE , x[4] , 99 ) ( VARIABLE , x[4] , 99 ) 
  ( CONSTANT , 3.3333333333333332e-04 , 99 ) 
  ( CONSTANT , 0.000212585 , 78 ) ( CONSTANT , -0.5 , 50 )
End
