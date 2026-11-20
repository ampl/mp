var x{1..5} integer >= -1 <= 1;

maximize o:
   6*x[1] + 4*x[2] - 2*x[3] + 2*x[4] - 4*x[5]
   - 2*(if abs(x[1] - x[5]) == 0 then 1 else 0);

s.t. c1: (x[4] != 0) <==> (x[4] != 0);
s.t. c2: (if x[4] != 0 then x[5] else 0) + (if x[4] != 0 then x[3] else 0) < -1;
s.t. c3: (abs(x[1] - x[5]) == 0) ==> (x[4] != 0) else (x[4] != 0);
