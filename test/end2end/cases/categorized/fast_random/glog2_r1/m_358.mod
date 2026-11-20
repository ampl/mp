var x{1..7} integer >= 0 <= 1;

maximize o: 6*x[1] + 6*x[2] - 4*x[3] + 0*x[4] + 6*x[5] - 4*x[6] + x[7];

s.t. c1: atleast 2 {i in 1..7} (x[i] >= 0);
s.t. c2: (x[2] == 0) and (x[4] > 1) ==> (x[7] >= 0);
