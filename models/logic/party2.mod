
param B > 0, integer;
set BOATS := 1 .. B;

param capacity {BOATS} integer >= 0;
param crew {BOATS} integer > 0;
param guest_cap {i in BOATS} := capacity[i] less crew[i];

param T > 0, integer;
set TIMES := 1..T;

set VISIT := {i in BOATS, j in BOATS: i <> j};
set MEET := {j in BOATS, jj in BOATS: j < jj};

       # (i,j) in VISIT if i can be visited by j
       # (j,jj) in MEET if j and jj can meet -- each pair only once

var isH {BOATS} binary;                   # i is a host boat

var H {BOATS,TIMES} integer, >= 1, <= B;  # boat visited by j at t

minimize TotalHosts: sum {i in BOATS} isH[i];

       # minimize total host boats

subj to VisitHosts {i in BOATS}:
   isH[i] = 0 ==> atmost 0 {j in BOATS, t in TIMES} (H[j,t] = i);

       # a visitor may not serve as a host: H[j,t] in VISITORS

subj to VisitOnce {j in BOATS}:
   isH[j] = 0 ==> alldiff {t in TIMES} H[j,t];
subj to HostNever {j in BOATS}:
   isH[j] = 1 ==> forall {t in TIMES} H[j,t] = j;

       # a non-host crew may visit a host at most once,
       # while a host crew does not visit at all

subj to Cap {i in BOATS, t in TIMES}: (isH[i] = 1) ==> 
   sum {(i,j) in VISIT} (if H[j,t] = i then crew[j]) <= guest_cap[i];

       # boats may not have more visitors than they can handle

subj to MeetDefn {(j,jj) in MEET}: (isH[j] = 0 and isH[jj] = 0) ==>
   atmost 1 {t in TIMES} (H[j,t] = H[jj,t]);

       # two crews may meet at most once

set MUST_BE_HOST within BOATS;
subj to MustBeHost {i in MUST_BE_HOST}: isH[i] = 1;

       # some boats are designated host boats

set MUST_BE_GUEST within BOATS;
subj to MustBeGuest {i in MUST_BE_GUEST}: isH[i] = 0;

       # some boats (the virtual boats) are designated guest boats

param mincrew := min {j in BOATS} crew[j];
subj to NeverHost {i in BOATS: guest_cap[i] < mincrew}: isH[i] = 0;

       # boats with very limited guest capacity can never be hosts
