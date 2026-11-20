var x{1..6} integer >= 0 <= 3;

maximize o: 0*x[1] + 6*x[2] + 3*x[3] - 5*x[4] - 2*x[5] + 0*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 3);
s.t. c2: abs(x[4] - x[3]) > 1;
