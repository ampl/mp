var x{1..5} integer >= 0 <= 4;

maximize o: 4*x[1] + 2*x[2] + 6*x[3] - 4*x[4] - 2*x[5];

s.t. c1: atleast 2 {i in 1..5} (x[i] >= 3);
s.t. c2: (x[5] < 0) ==> (x[5] != 4);
