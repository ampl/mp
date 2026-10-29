## Test multiple Q(C)P defvar substitution

param NX default 53;
param NDX default 10;
param NDXX default 25;

param lbx default -25.6;
param ubx default 13.4;
param lba default -3;
param uba default 6.3;
param lbb default -30;
param ubb default 50;

param p_inc_a default 0.6;
param p_inc_DX2 default 0.093;
param p_inc_DXX default 0.03;

param NSDXX_LB default 5;
param NSDXX_UB default 5;
param SDXX_LB default -500;
param SDXX_UB default 1500;

var x {1..NX} >=lbx <=ubx;

var dx {1..NDX}
   = sum {i in 1..NX}
      (if Uniform01() < p_inc_a
         then Uniform(lba, uba) * x[i])
     + Uniform(lbb, ubb);

var dxx {1..NDXX}
   = sum {i in 1..NDX, j in 1..NDX}
      if Uniform01() < p_inc_DX2
         then Uniform(0.5, 3.74) * dx[i] * dx[j]
     + Uniform(-13, 70.34);   ## This and/or the factor
          ## are essential for Gurobi 12.0.0 to have
          ## advantage with cvt:dvelim=2

minimize Obj:
   sum {i in 1..NDXX}
      if Uniform01() < p_inc_DXX
         then dxx[i];

s.t. Impose_SDXX_LB {k in 1..NSDXX_LB}:
   sum {i in 1..NDXX}
      if Uniform01() < p_inc_DXX
         then dxx[i]
   >= SDXX_LB;

s.t. Impose_SDXX_UB {k in 1..NSDXX_UB}:
   sum {i in 1..NDXX}
      if Uniform01() < p_inc_DXX
         then dxx[i]
   <= SDXX_UB;