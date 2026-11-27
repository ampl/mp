var x{1..3} integer >= 0 <= 3;

minimize o: -5*x[1] + x[2] + 0*x[3];

s.t. c1: atmost 3 {i in 1..3} (x[i] >= 2);
s.t. c2: x[2] == 1;
