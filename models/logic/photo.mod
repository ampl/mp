
# A group of people wants to take a group photo. Each person can give
# preferences next to whom he or she wants to be placed on the
# photo. The problem to be solved is to find a placement that
# satisfies as many preferences as possible.

# Adapted from
# www.g12.csse.unimelb.edu.au/minizinc/downloads/examples-latest/photo.mzn

param nPeople integer > 0;
set PREFS within {i1 in 1..nPeople, i2 in 1..nPeople: i1 <> i2};

var Sat {PREFS} binary;
var Pos {1..nPeople} integer >= 1, <= nPeople;

maximize NumSat: sum {(i1,i2) in PREFS} Sat[i1,i2];

subject to OnePersonPerPosition:
   alldiff {i in 1..nPeople} Pos[i];

subject to SatDefn {(i1,i2) in PREFS}:
   Sat[i1,i2] = 1 ==> Pos[i1]-Pos[i2] = 1 or Pos[i2]-Pos[i1] = 1;

subject to SymmBreaking:
   Pos[1] < Pos[2];
