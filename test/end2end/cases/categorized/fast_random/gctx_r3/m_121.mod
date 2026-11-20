var x{1..5} integer >= 0 <= 3;

maximize o:
   0*x[1] + 0*x[2] + 6*x[3] - 3*x[4] + 2*x[5] + (if x[3] >= 2 then 1 else 0);

s.t. c1: not (min(x[2], x[2]) < 1);
s.t. c2: not (min(x[2], x[2]) < 1);
s.t. c3: (max(x[3], x[5]) < 0) ==> (min(x[2], x[2]) < 1);
s.t. c4: not ((min(x[2], x[2]) < 1) and (min(x[2], x[2]) < 1));
