var x{1..4} integer >= -2 <= 1;

maximize o: -x[1] + 0*x[2] - 2*x[3] + 4*x[4];

s.t. c1: (x[2] >= -2) and (x[2] < 1) ==> (x[4] >= -2);
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 0);
