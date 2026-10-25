
# -------------------------------------------------------------
# CAR SEQUENCING (CSPLib prob001)
# carseq.mod: order the cars on an assembly line so that, for each
# option o, no more than p[o] cars in any window of q[o]
# consecutive slots require it (capacity p[o]/q[o]), while all the
# demands per class are met.
# Uses "numberof" (demands) and "count" (capacity windows).
# -------------------------------------------------------------

set CLASSES;
set OPTIONS;

param demand {CLASSES} integer >= 0;
param req {CLASSES, OPTIONS} binary;
param p {OPTIONS} integer > 0;
param q {OPTIONS} integer > 0;

param N := sum {c in CLASSES} demand[c];

set SLOTS := 1..N;

var Slot {SLOTS} integer in CLASSES;

subj to Demand {c in CLASSES}:
   numberof c in ({i in SLOTS} Slot[i]) = demand[c];

subj to Capacity {o in OPTIONS, s in 1..N-q[o]+1}:
   count {i in s..s+q[o]-1, c in CLASSES: req[c,o] = 1} (Slot[i] = c)
      <= p[o];

# Objective: get the cars of class 6 out as early as possible
minimize EarlyClass6:
   sum {i in SLOTS} i * (Slot[i] = 6);
