var x{1..4} integer >= 0 <= 3;

maximize o:
   5*x[1] + 6*x[2] + 0*x[3] - 5*x[4] + 0*(if x[4] != x[1] then 1 else 0);

s.t. c1:
   (if x[4] != x[1] then x[4] else 0) + (if x[4] != x[1] then x[2] else 0) >= 0;
s.t. c2: (max(x[3], x[3]) >= 3) ==> (max(x[3], x[3]) >= 3);
