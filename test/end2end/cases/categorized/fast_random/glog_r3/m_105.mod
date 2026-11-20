var x{1..3} integer >= 0 <= 4;

minimize o: 2*x[1] + 4*x[2] - 2*x[3];

s.t. c1: (x[2] <= 3) and (x[3] != 1) ==> (x[1] >= 3);
s.t. c2: (if x[1] == 3 then x[3] else x[1]) > 3;
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 0);
