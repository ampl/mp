
param B > 0, integer;
set BOATS := 1 .. B;
set HOSTS within BOATS;
set VISITORS := BOATS diff HOSTS;

param capacity {BOATS} integer >= 0;
param crew {BOATS} integer > 0;
param guest_cap {i in BOATS} := capacity[i] less crew[i];

param T > 0, integer;
set TIMES := 1..T;

set MEET := {j in VISITORS, jj in VISITORS: j < jj};

       # (j,jj) in MEET if j and jj can meet -- each pair only once

var H {VISITORS,TIMES} integer, >= 1, <= B;  # boat visited by j at t


VisitHosts {j in VISITORS, jj in VISITORS, t in TIMES}:
   H[j,t] != jj;

       # a visitor may not serve as a host: H[j,t] in VISITORS

subj to VisitOnce {j in VISITORS}:
   alldiff {t in TIMES} H[j,t];

       # a crew may visit each host at most once

subj to Cap {i in HOSTS, t in TIMES}: 
   sum {j in VISITORS} (if H[j,t] = i then crew[j]) <= guest_cap[i];

       # boats may not have more visitors than they can handle

subj to MeetDefn {(j,jj) in MEET}:
   atmost 1 {t in TIMES} (H[j,t] = H[jj,t]);

       # two crews may meet at most once
