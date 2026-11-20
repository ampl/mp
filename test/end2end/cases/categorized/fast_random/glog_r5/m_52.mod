var x{1..3} integer >= 0 <= 3;

minimize o: 5*x[1] - 5*x[2] - x[3];

s.t. c1: (x[1] != 0) ==> (x[1] < 2);
s.t. c2: atmost 1 {i in 1..3} (x[i] >= 1);
s.t. c3: not ((x[3] >= 2) and (x[2] >= 2));
