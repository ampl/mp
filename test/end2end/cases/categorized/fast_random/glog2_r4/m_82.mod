var x{1..4} integer >= 0 <= 3;

minimize o: -5*x[1] + 2*x[2] + 2*x[3] - 4*x[4];

s.t. c1: (x[3] >= 0) or (x[1] < 0) or (x[3] == 3);
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 3);
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 1);
