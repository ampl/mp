
set STU ordered;
param car {STU} >= 0, <= 1;

param ngroup integer >= 0;
set GRP := 1..ngroup;

set MEM {GRP} ordered by STU;
   check {g1 in GRP, g2 in g1+1..ngroup}: 
      card (MEM[g1] inter MEM[g2]) = 0;

set SAMEGRP := union {g in GRP}
   {s1 in MEM[g], s2 in MEM[g]: ord(s1) < ord(s2)};

# ----------------------------

set PRJ;
param cars_needed {PRJ} integer >= 0;
param min_team {PRJ} integer >= 0;
param max_team {p in PRJ} integer >= min_team[p];

param rank {STU,PRJ} integer >= 0, <= card {PRJ};

   check {(s1,s2) in SAMEGRP, p in PRJ}: rank[s1,p] = rank[s2,p];

# ----------------------------

var Assign {STU,PRJ} binary;

# ----------------------------

set GROUPED := union {g in GRP} MEM[g];

param group_weight >= 1;

# ----------------------------

minimize Total_Rank:
   sum {s in STU, p in PRJ} 
      (if s in GROUPED then 1 else group_weight) * rank[s,p] * Assign[s,p];

minimize Total_Rank_Squared:
   sum {s in STU, p in PRJ}  
      (if s in GROUPED then 1 else group_weight) * rank[s,p]^2 * Assign[s,p];

subject to Assign_Students {s in STU}:
   sum {p in PRJ} Assign[s,p] = 1;

subject to Assign_Projects {p in PRJ}:
   min_team[p] <= sum {s in STU} Assign[s,p] <= max_team[p];

subject to Enough_Cars {p in PRJ}:
   sum {s in STU} car[s] * Assign[s,p] >= cars_needed[p];

subject to Preserve_Groups {(s1,s2) in SAMEGRP, p in PRJ}:
   Assign[s1,p] = Assign[s2,p];

# ----------------------------

param cutoff >= 1, <= card {PRJ};

subject to Not_Too_Bad1 {s in STU, p in PRJ: rank[s,p] > cutoff}:
   Assign[s,p] = 0;

# ----------------------------

set PRJ_PREF within {STU,PRJ};

subject to Project_Preference {(s,p) in PRJ_PREF}: 
   Assign[s,p] = 1;

# ----------------------------

#set NO_RESPONSE := {s in STU: sum {p in PRJ} rank[s,p] = card(PRJ)};
#
#subject to Compensation {s0 in NO_RESPONSE, p in PRJ}:
#   sum {s in STU} rank[s,p] * Assign[s,p] 
#      <= 7 + card(PRJ) * (max {q in PRJ} max_team[q]) * (1 - Assign[s0,p]);
