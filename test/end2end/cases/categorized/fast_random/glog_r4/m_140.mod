var x{1..5} integer >= 0 <= 4;

maximize o: 3*x[1] - x[2] - 2*x[3] + 0*x[4] + 3*x[5];

s.t. c1: (x[4] < 4) and (x[4] == 0) ==> (x[2] >= 4);
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 0);
