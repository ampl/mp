var x{1..4} integer >= 0 <= 2;

minimize o:
   -x[1] + 6*x[2] + 6*x[3] + x[4] + 3*(if x[1] + x[3] > 0 then 1 else 0);

s.t. c1: (x[1] + x[3] > 0) and ((min(x[3], x[1]) < 0) or not (x[2] > 0));
s.t. c2: (x[4] <= 0) <==> (x[1] + x[3] > 0);
