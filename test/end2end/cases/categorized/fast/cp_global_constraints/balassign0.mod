
# ------------------------------------------------
# ANTI-ASSIGNMENT PROBLEM:
# Assign people to groups
# so that groups are as heterogeneous as possible
# ------------------------------------------------
# Version 2: Minimize "sum of deviations"
# ------------------------------------------------

# To make this problem harder,
# decrease sample and/or increase numberGrps.

set ALL_PEOPLE ordered;

param sample integer > 0;
param selection integer >= 0, < sample;
set PEOPLE := {i in ALL_PEOPLE: ord(i) mod sample = selection};

set CATEG;
param type {ALL_PEOPLE,CATEG} symbolic;
param typeWt {CATEG} >= 0;

param numberGrps integer > 0;

set TYPES {k in CATEG} := setof {i in PEOPLE} type[i,k];

var Assign {i in PEOPLE, j in 1..numberGrps} binary;

var MinInGrp <= floor (card(PEOPLE)/numberGrps) integer;
var MaxInGrp >= ceil (card(PEOPLE)/numberGrps) integer;

var MinType {k in CATEG, t in TYPES[k]} integer
   <= floor (card {i in PEOPLE: type[i,k] = t} / numberGrps);

var MaxType {k in CATEG, t in TYPES[k]} integer
   >= ceil (card {i in PEOPLE: type[i,k] = t} / numberGrps);

minimize Variation:  (MaxInGrp - MinInGrp) +
   sum {k in CATEG, t in TYPES[k]} 
      typeWt[k] * (MaxType[k,t] - MinType[k,t]);

subj to AssignAll {i in PEOPLE}:
   sum {j in 1..numberGrps} Assign[i,j] = 1;

subj to MinInGrpDefn {j in 1..numberGrps}:  
   MinInGrp <= sum {i in PEOPLE} Assign[i,j];

subj to MaxInGrpDefn {j in 1..numberGrps}:  
   MaxInGrp >= sum {i in PEOPLE} Assign[i,j];

subj to MinTypeDefn {j in 1..numberGrps, k in CATEG, t in TYPES[k]}:
   MinType[k,t] <= sum {i in PEOPLE: type[i,k] = t} Assign[i,j];

subj to MaxTypeDefn {j in 1..numberGrps, k in CATEG, t in TYPES[k]}:
   MaxType[k,t] >= sum {i in PEOPLE: type[i,k] = t} Assign[i,j];
