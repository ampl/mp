######################################################
## auxint_01.mod: implied-integer auxiliary variables.
## The nested if-then-else results (startup charge) and
## the count result are integer-valued given the binaries;
## with cvt:pre:auxint=0 (default) MP submits them as
## continuous and keeps their integrality internally.
######################################################

set G := 1..3;
set T := 1..4;

param cost {g in G} := 10*g;
param hot {g in G} := 50*g;
param cold {g in G} := 90*g;
param dem {t in T} := 120 + 20*t;

var u {G, T} binary;
var p {g in G, t in T} >= 0, <= 100;

var sc {g in G, t in T} =
    if t > 1 then
      (if u[g,t] = 1 and u[g,t-1] = 0 then
         (if t > 2 then (if u[g,t-2] = 1 then hot[g] else cold[g])
          else cold[g])
       else 0)
    else 0;

minimize Cost: sum {g in G, t in T} (cost[g]*p[g,t] + sc[g,t]);

s.t. Demand {t in T}: sum {g in G} p[g,t] >= dem[t];
s.t. Cap {g in G, t in T}: p[g,t] <= 100*u[g,t];
s.t. NStartups:
    count {g in G, t in T: t > 1} (u[g,t] = 1 and u[g,t-1] = 0) <= 2;
