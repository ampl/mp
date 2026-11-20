var x{1..4} integer >= -3 <= 1;

maximize o: 0*x[1] - 4*x[2] - x[3] + 0*x[4];

s.t. c1: atleast 1 {i in 1..4} (x[i] >= 1);
s.t. c2: (x[1] <= 1) and (x[3] > -1) ==> (x[1] >= 1);
