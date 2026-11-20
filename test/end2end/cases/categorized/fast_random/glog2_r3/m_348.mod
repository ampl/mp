var x{1..6} integer >= 0 <= 1;

maximize o: -4*x[1] + 4*x[2] + 5*x[3] + 0*x[4] - 4*x[5] - 3*x[6];

s.t. c1: atleast 1 {i in 1..6} (x[i] >= 0);
s.t. c2: ((if x[5] > 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[2] - x[5]) != 1;
s.t. c4: (x[4] > 0) and (x[5] > 0) ==> (x[1] >= 0);
