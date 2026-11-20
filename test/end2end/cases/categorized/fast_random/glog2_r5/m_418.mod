var x{1..7} integer >= 0 <= 1;

maximize o: -5*x[1] + 2*x[2] + 3*x[3] + 5*x[4] - 2*x[5] + 4*x[6] + 2*x[7];

s.t. c1: (x[2] <= 1) <==> (x[3] + x[2] != 1);
s.t. c2: atleast 3 {i in 1..7} (x[i] >= 1);
