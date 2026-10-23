\Nonlinear LP file
\Written from user model (presolvestate=0x680061)
\Names escaped using: @

\Problem name:                                                                  
\FICO Xpress v9.8.1, Hyper, written 16:06:44, May 18, 2026

Minimize
__OBJ___: - z + x - y + 0 LC1_flatv3_ + 0 LC1_flatv5_ + 0 LC1_flatv7_ + 0 LC1_flatv9_
 + 0 LC1_flatv10_ + 0 LC1_flatv11_ + 0 LC1_flatv12_ + 0 LC1_flatv13_

Subject To
R1: LC1_flatv2_=1 ->  3 z - x = -19 
R2: LC1_flatv4_=1 ->  2 z - 2 y = -27.99999 
R3: LC1_flatv6_=1 ->  x + 5 y = 65 
R4: LC1_flatv12_=1 ->  z = 3.8 
R5: LC1_flatv1_=1 ->  x = 0 
R6: LC1__condlineq4__auxv=1 ->  x + 5 y <= 64 
R7: LC1_flatv10_=1 ->  z + y <= 2.4 
R8: LC1__condlineq2__auxv=1 ->  3 z - x >= -18.9999 
R9: LC1__condlineq3__auxv=1 ->  2 z - 2 y >= -27.99989 
R10: LC1__condlineq1__auxv=1 ->  x >= 1 
R11: LC1_flatv2_ + LC1__condlineq2__auxv >= 1 
R12: LC1_flatv4_ + LC1__condlineq3__auxv >= 1 
R13: LC1_flatv6_ + LC1__condlineq4__auxv >= 1 
R14: z + x + y >= 15 
R15: LC1_flatv1_ + LC1__condlineq1__auxv >= 1 
R16: - LC1_flatv8_ + { ( LC1_flatv6_ * -1 + 1 ) * -1 + 1 } <= 0 

General Constraints
GC1: LC1_flatv13_ = OR ( LC1_flatv1_ , LC1_flatv2_ , LC1_flatv4_ , 
         LC1_flatv8_ , LC1_flatv10_ , LC1_flatv12_ )

Bounds
-3 <= z <= 7 
x <= 10 
y <= 11 
LC1_flatv3_ = 0 
LC1_flatv5_ = 0 
LC1_flatv7_ = 0 
LC1_flatv9_ = 0 
LC1_flatv11_ = 0 
LC1_flatv13_ = 1 

Generals
x y 

Binaries
LC1_flatv1_ LC1_flatv2_ LC1_flatv4_ LC1_flatv6_ LC1_flatv8_ LC1_flatv10_ 
LC1_flatv12_ LC1__condlineq2__auxv LC1__condlineq3__auxv LC1__condlineq4__auxv 
LC1__condlineq1__auxv 
 
NLP data
  LC1_flatv8_ determined_by            R16  0 0

End
