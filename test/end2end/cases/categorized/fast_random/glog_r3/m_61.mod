var x{1..4} integer >= 0 <= 2;

minimize o: 4*x[1] - 4*x[2] + 0*x[3] + x[4];

s.t. c1: abs(x[2] - x[2]) >= 0;
s.t. c2: (x[3] == 1) or (x[3] < 1) or (x[3] != 1);
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 2);
s.t. c4: (x[1] < 0) ==> (x[2] != 1) else (x[4] > 1);
