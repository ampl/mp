var x{1..6} integer >= -2 <= 1;

minimize o: x[1] - 5*x[2] + 2*x[3] + x[4] + 4*x[5] + 0*x[6];

s.t. c1: atmost 2 {i in 1..6} (x[i] >= 1);
s.t. c2: abs(x[2] - x[4]) != 3;
s.t. c3: atmost 2 {i in 1..6} (x[i] >= 0);
s.t. c4: abs(x[4] - x[6]) == 3;
