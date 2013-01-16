
### QUADRATIC ASSIGNMENT PROBLEM FROM GOLDMAN SACHS ###

set ALL_PEOPLE ordered;

param top {ALL_PEOPLE} symbolic;
param loc {ALL_PEOPLE} symbolic;
param bot {ALL_PEOPLE} symbolic;
param que {ALL_PEOPLE} symbolic;

param sample integer > 0;
set PEOPLE := {i in ALL_PEOPLE: ord(i) mod sample = 0} ordered by ALL_PEOPLE;

param topWt >= 0;
param locWt >= 0;
param botWt >= 0;
param queWt >= 0;

param numberGrps integer > 0;
param minInGrp := floor (card(PEOPLE)/numberGrps);
param nMinInGrp := numberGrps - card{PEOPLE} mod numberGrps;

param maxOverlap {i in PEOPLE} = minInGrp;

set TYPES = setof {i in PEOPLE} (top[i],loc[i],bot[i],que[i]);
set TYPEpeople {(t1,t2,t3,t4) in TYPES} =
   {i in PEOPLE: top[i]=t1 and loc[i]=t2 and bot[i]=t3 and que[i]=t4} ordered by PEOPLE;

var Assign {PEOPLE,1..numberGrps} binary;
var Overlap {PEOPLE} >= 0, <= 4 * minInGrp, integer;

minimize TotalOverlap:
   sum {i in PEOPLE} Overlap[i];

subj to OverlapDefn {i in PEOPLE, j in 1..numberGrps}:
   Overlap[i] >= 
      topWt * 
        sum {i2 in PEOPLE diff {i}: top[i2] = top[i]} Assign[i2,j] +
      locWt * 
        sum {i2 in PEOPLE diff {i}: loc[i2] = loc[i]} Assign[i2,j] +
      botWt * 
        sum {i2 in PEOPLE diff {i}: bot[i2] = bot[i]} Assign[i2,j] +
      queWt * 
        sum {i2 in PEOPLE diff {i}: que[i2] = que[i]} Assign[i2,j]
    - maxOverlap[i] * (1 - Assign[i,j]);

subj to AssignAll {i in PEOPLE}:
   sum {j in 1..numberGrps} Assign[i,j] = 1;

subj to GroupSizeMin {j in 1..nMinInGrp}:
   sum {i in PEOPLE} Assign[i,j] = minInGrp;
subj to GroupSizeMax {j in nMinInGrp+1..numberGrps}:
   sum {i in PEOPLE} Assign[i,j] = minInGrp + 1;

