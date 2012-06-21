
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

param nprj;
set PRJ := 1..nprj;

param cars_needed {PRJ} integer >= 0;
param min_team {PRJ} integer >= 0;
param max_team {p in PRJ} integer >= min_team[p];

param rank {STU,PRJ} integer >= 0, <= card {PRJ};

   check {(s1,s2) in SAMEGRP, p in PRJ}: rank[s1,p] = rank[s2,p];

# ----------------------------

var Project {STU} integer >= 1, <= nprj;

# ----------------------------

set GROUPED := union {g in GRP} MEM[g];

param group_weight >= 1;

# ----------------------------

minimize Total_Rank:
   sum {s in STU}  
      (if s in GROUPED then 1 else group_weight) * rank[s,Project[s]];

subject to Project_Size {p in PRJ}:
   min_team[p] <= numberof p in ({s in STU} Project[s]) <= max_team[p];

subject to Enough_Cars {p in PRJ}:
   numberof p in ({s in STU: car[s] = 1} Project[s]) >= cars_needed[p];

subject to Preserve_Groups {(s1,s2) in SAMEGRP}:
   Project[s1] = Project[s2];

# ----------------------------

param cutoff >= 1, <= card {PRJ};

 subject to Not_Too_Bad1 {s in STU, p in PRJ: rank[s,p] > cutoff}:
    Project[s] != p;

# ----------------------------

# set PRJ_PREF within {STU,PRJ};

# subject to Project_Preference {(s,p) in PRJ_PREF}: 
#    Project[s] = p;
