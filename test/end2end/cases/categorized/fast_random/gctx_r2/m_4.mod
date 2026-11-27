var x{1..5} integer >= 0 <= 3;

maximize o:
   0*x[1] - x[2] - x[3] + 4*x[4] + 0*x[5]
   + 2*(if x[5] + x[4] >= 1 then 1 else 0);

s.t. c1: not ((x[5] + x[5] == 1) and (x[5] + x[5] == 1));
s.t. c2: (x[5] + x[5] == 1) ==> (x[5] + x[4] >= 1) else (x[5] + x[4] >= 1);
