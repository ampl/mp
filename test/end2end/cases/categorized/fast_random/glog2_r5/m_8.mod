var x{1..7} integer >= 0 <= 1;

maximize o: 0*x[1] + 6*x[2] - 3*x[3] - 2*x[4] + 6*x[5] - 3*x[6] - 3*x[7];

s.t. c1: atmost 1 {i in 1..7} (x[i] >= 1);
s.t. c2: abs(x[6] - x[1]) > 1;
