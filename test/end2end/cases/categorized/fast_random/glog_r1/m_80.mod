var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] + 0*x[2] - 5*x[3];

s.t. c1: atmost 1 {i in 1..3} (x[i] >= 0);
s.t. c2: not ((x[2] != 2) and (x[3] >= 2));
s.t. c3: abs(x[2] - x[1]) < 0;
