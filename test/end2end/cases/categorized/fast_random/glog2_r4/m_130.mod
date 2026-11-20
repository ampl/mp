var x{1..3} integer >= -2 <= 1;

minimize o: -5*x[1] - x[2] + 5*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 1);
s.t. c2: (x[2] < 1) <==> (x[1] + x[2] < -2);
s.t. c3: not ((x[2] != -1) and (x[3] < -1));
