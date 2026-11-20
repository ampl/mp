var x{1..4} integer >= -1 <= 1;

maximize o: 0*x[1] - 4*x[2] - x[3] + 0*x[4] - 2*(if x[4] >= x[1] then 1 else 0);

s.t. c1: not (x[4] + x[3] >= 0);
s.t. c2: not (x[4] + x[3] >= 0);
s.t. c3: (x[4] != 1) ==> (min(x[2], x[3]) < 1);
