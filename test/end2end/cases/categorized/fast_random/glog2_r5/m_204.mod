var x{1..4} integer >= -3 <= 1;

maximize o: -x[1] + 6*x[2] + 4*x[3] + 0*x[4];

s.t. c1: (x[4] > -1) and (x[2] > -2) ==> (x[1] >= -1);
s.t. c2: atleast 1 {i in 1..4} (x[i] >= -1);
s.t. c3: ((if x[1] <= -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
