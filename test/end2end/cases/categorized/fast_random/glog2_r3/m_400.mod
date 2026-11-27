var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] + 2*x[2] - x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 2);
s.t. c2: max(x[1], x[3]) - min(x[3], x[2]) < 0;
