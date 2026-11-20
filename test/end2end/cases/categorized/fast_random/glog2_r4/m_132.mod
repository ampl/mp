var x{1..4} integer >= 0 <= 2;

minimize o: 2*x[1] - 4*x[2] - 2*x[3] + 4*x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= 0);
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: abs(x[1] - x[4]) >= 1;
