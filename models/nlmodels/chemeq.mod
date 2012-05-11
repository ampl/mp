# chemeq.mod	OLR2-MY-38-50
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Chemical equilibrium problem

# Ref.: D. M. Himmelblau, Applied Nonlinear Programming,
# McGraw-Hill Book Company, New York, 1972.  Problem 6.

# Number of variables:  38 (40 before presolve)
# Number of constraints:  50 (56 before presolve)
# Objective partially separable
# Objective convex
# Linear constraints

# Since constraints 14 and 16 imply x[1,5]=x[2,5]=x[3,5]=x[1,7]=x[2,7]=0
# these variables were removed from the original formulation of the problem

set J := 1 .. 18;
set K := 1 .. 7;
set PAIRS in {J,K};
set I := 1 .. 16;
param b {I} >= 0, default 0;
param c {PAIRS};
param E {PAIRS,I} integer;
param xlb >= 0, default 0;

var x {PAIRS} >= xlb, default 0.1;

minimize energy:
         sum {(j,k) in PAIRS} x[j,k] * (c[j,k] +
         log(x[j,k] / sum {m in J: (m,k) in PAIRS} x[m,k])) ;

s.t. h {i in I}:
     sum {(j,k) in PAIRS} E[j,k,i] * x[j,k] - b[i] = 0;

data;

param xlb := 1.5e-6;
set PAIRS:    1      2     3     4     5     6     7  :=
       1      +      +     +     +     -     +     -
       2      +      +     +     +     -     +     -
       3      +      +     +     +     -     -     -
       4      +      +     +     -     -     -     -
       5      -      +     +     -     -     -     -
       6      -      +     +     -     -     -     -
       7      -      +     +     -     -     -     -
       8      -      +     +     -     -     -     -
       9      -      +     +     -     -     -     -
      10      -      +     +     -     -     -     -
      11      -      +     +     -     -     -     -
      12      -      +     +     -     -     -     -
      13      -      +     +     -     -     -     -
      14      -      -     +     -     -     -     -
      15      -      -     +     -     -     -     -
      16      -      -     +     -     -     -     -
      17      -      -     +     -     -     -     -
      18      -      -     +     -     -     -     - ;

param      b         :=
      1   0.6529581
      2   0.281941
      3   3.705233
      4  47.00022
      5  47.02972
      6   0.08005
      7   0.08813
      8   0.04829
      9   0.0155
     10   0.0211275
     11   0.0022725;

param   c:     1        2        3        4        5        6        7  :=
        1     0      -10.94    10.45   -15.639    .          0        .
        2    -7.69     0        0        0        .        11.959     .
        3   -11.52     0       -0.5     21.81     .          .        .
        4   -36.60     0        0         .        .         .        .
        5      .       0        0         .        .         .        .
        6      .       0        0         .        .         .        .
        7      .       0        2.2435    .        .         .        .
        8      .       2.5966   0         .        .         .        .
        9      .     -39.39   -39.39      .        .         .        .
       10      .     -21.35   -21.49      .        .         .        .
       11      .     -32.84   -32.84      .        .         .        .
       12      .       6.26     6.12      .        .         .        .
       13      .       0        0         .        .         .        .
       14      .       .        0         .        .         .        .
       15      .       .       -1.9028    .        .         .        .
       16      .       .       -2.8889    .        .         .        .
       17      .       .       -3.3622    .        .         .        .
       18      .       .       -7.4854    .        .         .        .  ;


param           E   default 0   :=

         1 1 1      1
         2 1 2      1
         3 1 3      1
         4 1 4      1
         4 1 5      1
         1 2 1      1
         2 2 2      1
         3 2 3      1
         4 2 4      1
         4 2 12     1
         5 2 5      1
         5 2 12    -1
         6 2 6      1
         6 2 12    -1
         7 2 7      1
         7 2 12     1
         8 2 8      1
         8 2 12     1
         9 2 4      1
         9 2 5      1
        10 2 2      1
        10 2 5      1
        10 2 12    -1
        11 2 2      1
        11 2 4      1
        11 2 5      1
        12 2 2      1
        12 2 4     -1
        12 2 5      1
        12 2 12    -2
        13 2 9      1
        13 2 12    -1
         1 3 1      1
         2 3 2      1
         3 3 3      1
         4 3 4      1
         5 3 5      1
         6 3 6      1
         7 3 7      1
         8 3 8      1
         9 3 4      1
         9 3 5      1
        10 3 2      1
        10 3 5      1
        11 3 2      1
        11 3 4     -1
        11 3 5      1
        13 3 10     1
        14 3 11     1
        14 3 12    -4
        15 3 1      1
        15 3 11     1
        15 3 12    -3
        15 3 13    -1
        16 3 1      2
        16 3 11     1
        16 3 12    -2
        16 3 13    -2
        17 3 1      3
        17 3 11     1
        17 3 12    -1
        17 3 13    -3
        18 3 1      4
        18 3 11     1
        18 3 13    -4
         1 4 4      1
         1 4 13     1
         2 4 13     1
         3 4 4     -1
         3 4 13     1
         3 4 15    -4
         1 6 15     1
         2 6 2      1
         2 6 4     -1
         2 6 15     1  ;
