var x{1..3} integer >= -2 <= 1;

maximize o: -5*x[1] + 2*x[2] + x[3];

s.t. c1: abs(x[3] - x[3]) == 1;
s.t. c2: atmost 2 {i in 1..3} (x[i] >= 1);
